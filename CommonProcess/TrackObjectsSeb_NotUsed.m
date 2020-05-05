function Var = TrackObjectsSeb(Var, CallNum)

Debug = 0;      %Set Debug to 0 to prevent image display

if nargin == 1
    CallNum = 1;
end
%Get Track Object and Tracking distance from Var
TrackObj = 'Nucl' ; Var.Analysis.TrackObj{CallNum};
TrackIllum = 'CorrCFP';
MaxDistance = Var.Analysis.TrackMaxDistance{CallNum};


%Loop through all frames
for F = 1:2 %Var.Analysis.NumFrame
    %TrackImg = zeros(Var.Analysis.ImgSize);
    
    %Read Frame.mat file
    ImgFolder = fileparts(Var.Analysis.OutPath);
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(F,'%05d'),'.mat'])
    load(FileName)
    
    %Get X and Y centers from measurement
    CX{F} = Var.Measurements.(TrackObj).(TrackIllum).CenterX';
    CY{F} = Var.Measurements.(TrackObj).(TrackIllum).CenterY';
    ObjPixel{F} = Var.Measurements.(TrackObj).PixelList;
    SegLabel{F} = Var.Measurements.(TrackObj).SegLabel;
    Area{F} = Var.Measurements.(TrackObj).(TrackIllum).Area;
    Orientation{F} = Var.Measurements.(TrackObj).(TrackIllum).Orientation;
    MajAxis{F} = Var.Measurements.(TrackObj).(TrackIllum).MajorAxisLength;
    MinAxis{F} = Var.Measurements.(TrackObj).(TrackIllum).MinorAxisLength;
    
    
    for Obj = 1:length(ObjPixel{F})
        ObjImg = zeros(Var.Analysis.ImgSize);
        ObjImg(ObjPixel{F}{Obj}) = 1;
        ObjImg = bwmorph(ObjImg,'remove');
        [PerimPixel{F}{Obj}(:,1), PerimPixel{F}{Obj}(:,2)] = find(ObjImg == 1);
    end
        
    
    
end

lambda(1) = 5; % Weight for Delta_pos
lambda(2) = 0.1; % Weight for weighted_posx
lambda(3) = 2; % Weight for weighted_posy
lambda(4) = 0.1; % Weight for Delta_area
lambda(5) = 50;
lambda(6) = 10;
maxfact = 9;

Previous = 1;
Current = 2;


single_cell_indices1 = SegLabel{Previous};
single_cell_indices2 = SegLabel{Current};

S_si1 = size(single_cell_indices1)
S_si2 = size(single_cell_indices2)

centroid1  = [CX{Previous}, CY{Previous}];
centroid2  = [CX{Current}, CY{Current}];

center1  = [CX{Previous}, CY{Previous}];
center2  = [CX{Current}, CY{Current}];

S_C1 = size(center1)
S_c2 = size(center2)

slope1 =  Orientation{Previous};
slope2 =  Orientation{Current};

boundarea1 = Area{Previous};
boundarea2 = Area{Current};

longueur1 =  MajAxis{Previous};
longueur2 =  MajAxis{Current};

larg1 = MinAxis{Previous};
larg2 = MinAxis{Current};


linkage_Map_pos = zeros(length(CX{Previous}), length(CX{Current}));
linkage_Map_area = zeros(length(CX{Previous}), length(CX{Current}));
linkage_posx = zeros(length(CX{Previous}), length(CX{Current}));
linkage_posy = zeros(length(CX{Previous}), length(CX{Current}));
linkage_posxreal = zeros(length(CX{Previous}), length(CX{Current}));
linkage_area_ratio = zeros(length(CX{Previous}), length(CX{Current}));
orientation_Map_area = zeros(length(CX{Previous}), length(CX{Current}));
scoring_Map = zeros(length(CX{Previous}), length(CX{Current}));
Test_variable = zeros(length(CX{Previous}), length(CX{Current}));
linkage_move = zeros(length(CX{Previous}), length(CX{Current}));
%We will calculate a deltax and a deltay which will be directly related to
%the angle of the cells + a "normal" DeltaX

for k = 1:length(CX{Previous})
    for l = 1:length(CX{Current})
        marqueur = 0;
        
        delta_posxreal = CX{Current}(l)-CX{Previous}(k);
        delta_posyreal = CY{Current}(l)-CY{Previous}(k);
        weighted_posx = delta_posxreal*0.8;
        
        if delta_posxreal < -10
            weighted_posx = delta_posxreal * 2;
        end
        
        
        area = Area{Previous}(k)-Area{Current}(l);
        delta_area = abs((Area{Previous}(k)-Area{Current}(l)));
        area_ratio = delta_area/max([Area{Previous}(k)  Area{Current}(l)]);
        
        angle2 = Orientation{Current}(l);
        angle1 = Orientation{Previous}(k);
        if (angle2 > pi/2-0.4 && angle2 < pi/2+0.4)
            if delta_posxreal < -2
                weighted_posx = delta_posxreal*2;
            end
        end
        
        
        delta_orientation = Delta_Angle(angle1,angle2);
        
        delta_pos = sqrt((weighted_posx)^2 + (delta_posyreal)^2);
        if delta_pos > MaxDistance
            scoring_Map(k,l) = Inf;
        else
            linkage_Map_pos(k,l) = delta_pos;
            linkage_posxreal(k,l)=delta_posxreal;
            linkage_Map_area(k,l) = area;
            orientation_Map_area(k,l) = delta_orientation;
            linkage_area_ratio(k,l) = area_ratio;
            scoring_Map(k,l) = lambda(1) * delta_pos + lambda(2) *abs(weighted_posx) + lambda(3) *abs(delta_posyreal) + lambda(4)*delta_area + lambda(5)* delta_orientation;
        end
        
    end
end

assignin('base', 'scoring_Map', scoring_Map)
           C = zeros(length(single_cell_indices2),1);
            I = zeros(length(single_cell_indices2),1);
            C2 = zeros(length(single_cell_indices2),1);
            I2 = zeros(length(single_cell_indices2),1);
            division_marker = zeros(length(single_cell_indices2),1);
            clear to_work;
for m = 1:length(single_cell_indices2)
    
    [sorted_min, ix_sorted] = sort(scoring_Map(:,m));
    voisins = ix_sorted(sorted_min ~= Inf);
    sorted_min(1:length(voisins));
    compt=1;
    boundary2 = PerimPixel{Current}{single_cell_indices2(m)};
    boundary1 = PerimPixel{Previous}{single_cell_indices1(ix_sorted(compt))};
    [x1,y1] = poly2cw(boundary1(:,2), boundary1(:,1));
    [x2,y2] = poly2cw(boundary2(:,2), boundary2(:,1));
    
    [xb, yb] = polybool('intersection', x1, y1, x2, y2);
    if isempty(xb) == 0
        
        if sorted_min(compt)<65
            if linkage_area_ratio(ix_sorted(compt), m) < 0.4
                C(m) = scoring_Map(ix_sorted(compt),m);
                I(m) = ix_sorted(compt);
                to_work{m} = 0;
            else
                to_work{m} = voisins;
            end
        else
            to_work{m} = voisins;
        end
    else
        to_work{m} = voisins;
    end
end
clear f2un;
f2un = find(I==0);
clear pair;
clear f1un
compt = 1;
for k = 1:length(single_cell_indices1)
    if isempty(find(I == k)) == 1
        f1un(compt) = k;
        compt = compt+1;
    end
end
clear f1unmatched
clear f2unmatched
compt = 1;
compt2 = 1;
max_move = 0;
total_displacement = 0;
for k = 1:length(single_cell_indices1)
    if isempty(find(I == k)) == 1;
        f1un(compt) = k;
        compt = compt+1;
    elseif length(find(I == k)) > 1 || isempty(find(I2 == k)) ~=1
        
    else
        
        distance = sqrt(sum(diff([centroid1(single_cell_indices1(k),:) ; centroid2(single_cell_indices2(find(I ==k)),:)]).^2));
        dis(compt2) = distance;
        compt2 = compt2+1;
        if distance > max_move
            max_move = distance;
            pos = k;
        end
        if size(distance,2) ~= 1
            'cry'
            k
            break
        end
        total_displacement = total_displacement + distance;
    end
end



mean_move = total_displacement/(compt2-1);
if exist('f1un','var') == 1 && exist('f2un', 'var') == 1
    f1un = f1un';
    match_compt = 1;
    clear pair;
    
    f1untemp = f1un
    f2untemp = f2un
    for l = 1:size(f1un,1)
        cvois = 0;
        clear lisvois;
        bound1 = PerimPixel{Previous}{single_cell_indices1(f1un(l))};
        cent1 =  center1(single_cell_indices1(f1un(l)),:);
        if cent1(1)>10
            angle1 = slope1(single_cell_indices1(f1un(l)));
            area1 = boundarea1(single_cell_indices1(f1un(l)));
            long1 = longueur1(single_cell_indices1(f1un(l)));
            SZ_f2 = size(f2un,1)
            
            for l2 = 1:size(f2un,1);
                if f2untemp(l2,1) ~= 0
                    bound2 = PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                    assignin('base', 'PerimPixel', PerimPixel)
                    assignin('base', 'center2', center2)
                    assignin('base', 'single_cell_indices2', single_cell_indices2)
                    assignin('base', 'f2un', f2un)
                    assignin('base', 'l2', l2)
                    cent2 = center2(single_cell_indices2(f2un(l2)),:);
                    long2 = longueur2(single_cell_indices2(f2un(l2)));
                    if  sqrt(sum(diff([cent1;cent2]).^2)) < max([long1/1.8 long2/1.8 mean_move])
                        delta_posxreal = cent2(1)-cent1(1);
                        delta_posyreal = cent2(2)-cent1(2);
                        if  delta_posxreal > 1 && delta_posyreal > 1
                            moving = atan(delta_posyreal/delta_posxreal);
                            if moving < 0
                                moving = moving+pi;
                            end
                        else
                            moving = 0;
                        end
                        angle2 = slope2(single_cell_indices2(f2un(l2)));
                        
                        delta_orientation = Delta_Angle(angle1,angle2);
                        if long2 < 25
                            larg_lon_ratio = longueur2(single_cell_indices2(f2un(l2)))/larg2(single_cell_indices2(f2un(l2)));
                            if larg_lon_ratio < 1.5  || cent1(1) > 1200
                                maxangle = pi/2;
                            else
                                maxangle = 0.7;
                            end
                        else
                            maxangle = 0.7;
                        end
                        if delta_orientation < maxangle
                            cvois = cvois+1;
                            lisvois(cvois) = l2;
                            area2 = boundarea2(single_cell_indices2(f2un(l2)));
                            polyarea(bound2(:,2), bound2(:,1));
                        end
                    end
                end
            end
            if cvois == 1
                l2 = lisvois(1);
                long2 = longueur2(single_cell_indices2(f2un(l2)));
                area2 =  boundarea2(single_cell_indices2(f2un(l2)));
                angle2 = slope2(single_cell_indices2(f2un(l2)));
                delta_orientation = Delta_Angle(angle1,angle2);
                fla = 0;
                if max([long1 long2]) > 40
                    if abs(long2-long1)/(long2+long1) < 0.2 && delta_orientation < 0.4
                        fla =1;
                    end
                else
                    if abs(long2-long1)/(long2+long1) < 0.35 && delta_orientation < 0.4
                        fla =1;
                    end
                end
                if fla == 1
                    
                    bound2 =  PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                    conv2 = convhull(bound2(:,2), bound2(:,1));
                    conarea2 = polyarea(bound2(conv2,2), bound2(conv2,1));
                    conv1 = convhull(bound1(:,2), bound1(:,1));
                    conarea1 = polyarea(bound1(conv1,2), bound1(conv1,1));
                    if abs(area2-area1)/(area2+area1) < 0.35 || abs(conarea2-conarea1)/(conarea2+conarea1) < 0.33
                        pair(match_compt,:) = [l l2];
                        match_compt = match_compt+1;
                        f1untemp(l,:) = 0;
                        f2untemp(l2,:) = 0;
                    end
                    
                end
            elseif cvois == 2
                clear bounda;
                clear centa;
                clear anglea;
                clear longa;
                clear disX
                clear disY
                for k = 1:cvois
                    l2 = lisvois(k);
                    bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                    centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                    anglea(k) = slope2(single_cell_indices2(f2un(l2)));
                    longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                    disX(k) = centa{k}(1)-cent1(1);
                    disY(k) = centa{k}(2)-cent1(2);
                end
                delta_orientation = Delta_Angle(anglea(k),anglea(k-1));
                if (delta_orientation < 0.30 && (bounda(k)+bounda(k-1))/area1 > 0.72 && (bounda(k)+bounda(k-1))/area1 < 1.28) || ((sum(longa)-long1)/(sum(longa)+long1) < 0.1 && (bounda(k)+bounda(k-1))/area1 > 0.8250 && (bounda(k)+bounda(k-1))/area1 < 1.175) ||  ((sum(longa)-long1)/(sum(longa)+long1) < 0.1 && delta_orientation < 0.30 || ((bounda(k)+bounda(k-1))/area1 < 1.1 &&  (bounda(k)+bounda(k-1))/area1 > 0.9 && delta_orientation < 0.45))
                    meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                    if abs(meancent(1)- cent1(1)) < max([mean_move 15]) && abs(meancent(2) - cent1(2)) < max([mean_move 11])
                        for k = 1:cvois
                            l2 = lisvois(k);
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = 0;
                            f2untemp(l2,:) = 0;
                            bound2 =  PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                        end
                    end
                elseif delta_orientation < 0.35
                    bound1 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois(1)))};
                    bound2 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois(2)))};
                    [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound2(:,2)], [bound1(:,1); bound2(:,1)]);
                    first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                    second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                    if first < second
                        width = first;
                        longue = second;
                    else
                        width = second;
                        longue = first;
                    end
                    longueur = longue;
                    if abs((longueur-long1)/(longueur+long1)) < 0.1 && delta_orientation < 0.35
                        meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                        if abs(meancent(1)- cent1(1)) < max([mean_move 15]) && abs(meancent(2) - cent1(2)) < max([mean_move 11])
                            for k = 1:cvois
                                l2 = lisvois(k);
                                pair(match_compt,:) = [l l2];
                                match_compt = match_compt+1;
                                f1untemp(l,:) = 0;
                                f2untemp(l2,:) = 0;
                                bound2 =  PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                            end
                        end
                    end
                else
                    ar2 = abs(bounda(k)-area1)/(area1 + bounda(k));
                    ar1 = abs(bounda(k-1)-area1)/(area1 + bounda(k-1));
                    if abs(ar1-ar2) > 0.15
                        [minar minaridx] = min([ar1, ar2]);
                        k = minaridx;
                        delta_or = Delta_Angle(anglea(k), angle1);
                        l2 = lisvois(k);
                        delta_posxreal = centa{k}(1)-cent1(1);
                        delta_posyreal = centa{k}(2)-cent1(2);
                        if abs(delta_posxreal) > 1 || abs(delta_posyreal) > 1
                            moving(k) = atan(delta_posyreal/delta_posxreal);
                        else
                            moving(k) = 0;
                        end
                        if moving(k) < 0
                            moving(k) = moving(k)+pi;
                        end
                        Delta_move = Delta_Angle(moving(k), angle1);
                        if ((abs(longa(k)-long1)/(longa(k)+long1) < 0.2 && abs(bounda(k)-area1)/(bounda(k)+area1) < 0.15) || (abs(longa(k)-long1)/(longa(k)+long1) < 0.3 && abs(bounda(k)-area1)/(bounda(k)+area1) < 0.3 && delta_or <0.15)) &&  disX(k) < 8 && Delta_move < 0.9
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = [0];
                            f2untemp(l2,:) = [0];
                        end
                    end
                    
                end
            elseif cvois > 2
                cvois2 = 0;
                clear lisvois2;
                clear bounda;
                clear centa;
                clear anglea;
                clear longa;
                for k = 1:cvois
                    l2 = lisvois(k);
                    bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                    centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                    anglea(k) = slope2(single_cell_indices2(f2un(l2)));
                    longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                    delta_posxreal = centa{k}(1)-cent1(1);
                    delta_posyreal = centa{k}(2)-cent1(2);
                    if long1 < 40
                        if abs(delta_posxreal) > 1 || abs(delta_posyreal) > 1
                            moving(k) = atan(delta_posyreal/delta_posxreal);
                            if moving(k) < 0
                                moving(k) = moving(k)+pi;
                            end
                        else moving(k) = 0;
                        end
                    else
                        if abs(delta_posxreal) > long1/30 || abs(delta_posyreal) > long1/30;
                            moving(k) = atan(delta_posyreal/delta_posxreal);
                            if moving(k) < 0
                                moving(k) = moving(k)+pi;
                            end
                        else moving(k) = 0;
                        end
                    end
                    delta_or = Delta_Angle(moving(k), angle1);
                    if Delta_Angle(moving(k), 0) < 0.4 || delta_or < 0.5 || (delta_posxreal < 3*mean_move && delta_or < 0.9 && abs(longa(k)-long1)/(longa(k)+long1) < 0.2)
                        cvois2 = cvois2+1;
                        lisvois2(cvois2) = l2;
                        
                    end
                end
                
                
                if cvois2 > 3
                    
                    cvois3 = 0;
                    clear lisvois3;
                    clear bounda;
                    clear centa;
                    clear anglea;
                    clear longa;
                    for k = 1:cvois2
                        l2 = lisvois2(k);
                        bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                        centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                        anglea(k) = slope2(single_cell_indices2(f2un(l2)));
                        longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                        delta_posxreal = centa{k}(1)-cent1(1);
                        delta_posyreal = centa{k}(2)-cent1(2);
                        if long1 < 40
                            if abs(delta_posxreal) > 1 || abs(delta_posyreal) > 1
                                moving(k) = atan(delta_posyreal/delta_posxreal);
                                if moving(k) < 0
                                    moving(k) = moving(k)+pi;
                                end
                            else moving(k) = 0;
                            end
                        else
                            if abs(delta_posxreal) > long1/30 || abs(delta_posyreal) > long1/30;
                                moving(k) = atan(delta_posyreal/delta_posxreal);
                                if moving(k) < 0
                                    moving(k) = moving(k)+pi;
                                end
                            else moving(k) = 0;
                            end
                        end
                        delta_or = Delta_Angle(moving(k), angle1);
                        if Delta_Angle(moving(k), 0) < 0.4 || delta_or < 0.3
                            cvois3 = cvois3+1;
                            lisvois3(cvois3) = l2;
                            
                        end
                    end
                    if cvois3 > 2
                        cvois4 = 0;
                        clear lisvois4;
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear delta_posxreal;
                        clear delta_posyreal;
                        clear moving
                        clear distance
                        for k = 1:cvois3
                            l2 = lisvois3(k);
                            bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                            centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                            anglea(k) =slope2(single_cell_indices2(f2un(l2)));
                            longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                            delta_posxreal(k) = centa{k}(1)-cent1(1);
                            delta_posyreal(k) = centa{k}(2)-cent1(2);
                            distance(k) = sqrt(sum(diff([centa{k};cent1]).^2));
                        end
                        bound1 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois3(1)))};
                        bound2 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois3(2)))};
                        bound3 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois3(3)))};
                        [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound2(:,2)], [bound1(:,1); bound2(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(1) = first;
                            longue(1) = second;
                        else
                            width(1) = second;
                            longue(1) = first;
                        end
                        
                        
                        [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound3(:,2)], [bound1(:,1); bound3(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(2) = first;
                            longueur(2) = second;
                        else
                            width(2) = second;
                            longue(2) = first;
                        end
                        
                        
                        [rectx, recty, aire, perim] =  minboundrect([bound3(:,2); bound2(:,2)], [bound3(:,1); bound2(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(3) = first;
                            longue(3) = second;
                        else
                            width(3) = second;
                            longue(3) = first;
                        end
                        [valmin minidx] = min(width);
                        k = minidx;
                        if k == 1
                            lisvois4 = [lisvois2(1) lisvois2(2)];
                        elseif k == 2
                            lisvois4 = [lisvois2(1) lisvois2(3)];
                        elseif k == 3
                            lisvois4 = [lisvois2(2) lisvois2(3)];
                        end
                        longueur = longue(k);
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear delta_posxreal;
                        clear delta_posyreal;
                        clear moving
                        clear distance
                        for k = 1:2
                            l2 = lisvois4(k);
                            bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                            centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                            anglea(k) =slope2(single_cell_indices2(f2un(l2)));
                            longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                            delta_posxreal(k) = centa{k}(1)-cent1(1);
                            delta_posyreal(k) = centa{k}(2)-cent1(2);
                            distance(k) = sqrt(sum(diff([centa{k};cent1]).^2));
                        end
                        
                        delta_or = Delta_Angle(anglea(k), anglea(k-1));
                        if (delta_or < 0.3 && ((bounda(k)+bounda(k-1))/area1 > 0.75 && (bounda(k)+bounda(k-1))/area1 < 1.25) || abs((longueur-long1))/(longueur+long1) < 0.2 ) || (delta_or < 0.8 && (bounda(k)+bounda(k-1))/area1 > 0.85 && (bounda(k)+bounda(k-1))/area1 < 1.15)
                            meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                            if abs(meancent(1)- cent1(1)) < max([abs(cos(angle1))*long1/2, max_move]) && abs(meancent(2) - cent1(2)) <  abs(sin(angle1))*long1/2+2
                                for k = 1:2
                                    l2 = lisvois4(k);
                                    pair(match_compt,:) = [l l2];
                                    match_compt = match_compt+1;
                                    f1untemp(l,:) = 0;
                                    f2untemp(l2,:) = 0;
                                end
                            end
                        end
                    elseif cvois3 == 2
                        clear bounda
                        clear centa
                        clear anglea
                        clear longa
                        clear moving;
                        for k = 1:cvois3
                            l2 = lisvois3(k);
                            bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                            centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                            anglea(k) = slope2(single_cell_indices2(f2un(l2)));
                            longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                        end
                        delta_orientation = Delta_Angle(anglea(k), anglea(k-1));
                        if ((bounda(k)+bounda(k-1))/area1 > 0.75 && (bounda(k)+bounda(k-1))/area1 < 1.25 && delta_orientation < 0.3 ) || ((bounda(k)+bounda(k-1))/area1 > 0.85 && (bounda(k)+bounda(k-1))/area1 < 1.15 && delta_orientation < 0.5 )
                            meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                            if abs(meancent(1)- cent1(1)) < abs(cos(angle1))*long1/2 && abs(meancent(2) - cent1(2)) <  max([abs(sin(angle1))*long1/2 2])
                                for k = 1:cvois3
                                    l2 = lisvois3(k);
                                    pair(match_compt,:) = [l l2];
                                    match_compt = match_compt+1;
                                    f1untemp(l,:) = 0;
                                    f2untemp(l2,:) = 0;
                                    bound2 =  PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                                end
                            end
                        end
                    elseif cvois3 == 1
                        l2 = lisvois3(1);
                        area2 = boundarea2(single_cell_indices2(f2un(l2)));
                        angle2 = slope2(single_cell_indices2(f2un(l2)));
                        long2 = longueur2(single_cell_indices2(f2un(l2)));
                        delta_orientation = Delta_Angle(angle1, angle2);
                        if abs(long2-long1)/(long2+long1) < 0.2 && delta_orientation < 0.2 && abs(area2-area1)/(area2+area1) < 0.1
                            bound2 = PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = 0;
                            f2untemp(l2,:) = 0;
                            
                            
                        end
                    end
                elseif cvois2 == 3
                    cvois3 = 0;
                    clear lisvois3;
                    clear bounda;
                    clear centa;
                    clear anglea;
                    clear longa;
                    clear delta_posxreal;
                    clear delta_posyreal;
                    clear moving
                    clear distance
                    for k = 1:cvois2
                        l2 = lisvois2(k);
                        bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                        centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                        anglea(k) =slope2(single_cell_indices2(f2un(l2)));
                        longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                        delta_posxreal(k) = centa{k}(1)-cent1(1);
                        delta_posyreal(k) = centa{k}(2)-cent1(2);
                        distance(k) = sqrt(sum(diff([centa{k};cent1]).^2));
                    end
                    bound1 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois2(1)))};
                    bound2 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois2(2)))};
                    bound3 = PerimPixel{Current}{single_cell_indices2(f2un(lisvois2(3)))};
                    
                    [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound2(:,2)], [bound1(:,1); bound2(:,1)]);
                    first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                    second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                    if first < second
                        width(1) = first;
                        longue(1) = second;
                    else
                        width(1) = second;
                        longue(1) = first;
                    end
                    
                    
                    [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound3(:,2)], [bound1(:,1); bound3(:,1)]);
                    first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                    second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                    if first < second
                        width(2) = first;
                        longueur(2) = second;
                    else
                        width(2) = second;
                        longue(2) = first;
                    end
                    
                    
                    [rectx, recty, aire, perim] =  minboundrect([bound3(:,2); bound2(:,2)], [bound3(:,1); bound2(:,1)]);
                    first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                    second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                    if first < second
                        width(3) = first;
                        longue(3) = second;
                    else
                        width(3) = second;
                        longue(3) = first;
                    end
                    [valmin minidx] = min(width);
                    k = minidx;
                    if k == 1
                        lisvois3 = [lisvois2(1) lisvois2(2)];
                    elseif k == 2
                        lisvois3 = [lisvois2(1) lisvois2(3)];
                    elseif k == 3
                        lisvois3 = [lisvois2(2) lisvois2(3)];
                    end
                    longueur = longue(k);
                    clear bounda;
                    clear centa;
                    clear anglea;
                    clear longa;
                    clear delta_posxreal;
                    clear delta_posyreal;
                    clear moving
                    clear distance
                    for k = 1:2
                        l2 = lisvois3(k);
                        bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                        centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                        anglea(k) =slope2(single_cell_indices2(f2un(l2)));
                        longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                        delta_posxreal(k) = centa{k}(1)-cent1(1);
                        delta_posyreal(k) = centa{k}(2)-cent1(2);
                        distance(k) = sqrt(sum(diff([centa{k};cent1]).^2));
                    end
                    
                    delta_or = Delta_Angle(anglea(k), anglea(k-1));
                    if (delta_or < 0.3 && ((bounda(k)+bounda(k-1))/area1 > 0.75 && (bounda(k)+bounda(k-1))/area1 < 1.25) || abs((longueur-long1))/(longueur+long1) < 0.2 ) || (delta_or < 0.8 && (bounda(k)+bounda(k-1))/area1 > 0.85 && (bounda(k)+bounda(k-1))/area1 < 1.15)
                        meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                        if abs(meancent(1)- cent1(1)) < max([abs(cos(angle1))*long1/2, max_move]) && abs(meancent(2) - cent1(2)) <  abs(sin(angle1))*long1/2+2
                            for k = 1:2
                                l2 = lisvois3(k);
                                pair(match_compt,:) = [l l2];
                                match_compt = match_compt+1;
                                f1untemp(l,:) = 0;
                                f2untemp(l2,:) = 0;
                            end
                        end
                    end
                elseif cvois2 == 2
                    clear bounda
                    clear centa
                    clear anglea
                    clear longa
                    clear moving;
                    for k = 1:cvois2
                        l2 = lisvois2(k);
                        bounda(k) = boundarea2(single_cell_indices2(f2un(l2)));
                        centa{k} = center2(single_cell_indices2(f2un(l2)),:);
                        anglea(k) = slope2(single_cell_indices2(f2un(l2)));
                        longa(k) = longueur2(single_cell_indices2(f2un(l2)));
                    end
                    delta_orientation = Delta_Angle(anglea(k), anglea(k-1));
                    if ((bounda(k)+bounda(k-1))/area1 > 0.75 && (bounda(k)+bounda(k-1))/area1 < 1.25 && delta_orientation < 0.3 ) || ((bounda(k)+bounda(k-1))/area1 > 0.85 && (bounda(k)+bounda(k-1))/area1 < 1.15 && delta_orientation < 0.5 ) || (sum(longa)/long1 < 1.15 && sum(longa)/long1 > 0.85 && delta_orientation < 0.5)
                        meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                        if abs(meancent(1)- cent1(1)) < abs(cos(angle1))*long1/2 && abs(meancent(2) - cent1(2)) <  max([abs(sin(angle1))*long1/2 2])
                            for k = 1:cvois2
                                l2 = lisvois2(k);
                                pair(match_compt,:) = [l l2];
                                match_compt = match_compt+1;
                                f1untemp(l,:) = 0;
                                f2untemp(l2,:) = 0;
                                bound2 =  PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                            end
                        end
                    end
                elseif cvois2 == 1
                    l2 = lisvois2(1);
                    area2 = boundarea2(single_cell_indices2(f2un(l2)));
                    angle2 = slope2(single_cell_indices2(f2un(l2)));
                    long2 = longueur2(single_cell_indices2(f2un(l2)));
                    delta_orientation = Delta_Angle(angle1, angle2);
                    if abs(long2-long1)/(long2+long1) < 0.2 &&delta_orientation < 0.4 && abs(area2-area1)/(area2+area1) < 0.1
                        bound2 = PerimPixel{Current}{single_cell_indices2(f2un(l2))};
                        pair(match_compt,:) = [l l2];
                        match_compt = match_compt+1;
                        f1untemp(l,:) = 0;
                        f2untemp(l2,:) = 0;
                    end
                end
                
                
                
            end
        end
    end
    for l2 = 1:size(f2un,1)
        %                 l2 = l2+1
        cvois = 0;
        clear lisvois;
        cvoisput = 0;
        clear lisvoisput;
        if f2untemp(l2,1) ~= 0
            bound2 = PerimPixel{Current}{single_cell_indices2(f2un(l2))};
            cent2 = center2(single_cell_indices2(f2un(l2)),:);
            if cent2(1)>30
                angle2 = slope2(single_cell_indices2(f2un(l2)));
                area2 = boundarea2(single_cell_indices2(f2un(l2)));
                long2 = longueur2(single_cell_indices2(f2un(l2)));
                for l = 1:size(f1un,1);
                    if f1untemp(l,1) ~= 0
                        bound1 = PerimPixel{Previous}{single_cell_indices1(f1un(l))};
                        cent1 =  center1(single_cell_indices1(f1un(l)),:);
                        long1 = longueur1(single_cell_indices1(f1un(l)));
                        if   sqrt(sum(diff([cent1;cent2]).^2)) < max([long1/2 long2/2 mean_move])
                            angle1 =slope1(single_cell_indices1(f1un(l)));
                            delta_posxreal = cent2(1)-cent1(1);
                            delta_posyreal = cent2(2)-cent1(2);
                            if  abs(delta_posxreal) > 1 && abs(delta_posyreal) > 1
                                moving = atan(delta_posyreal/delta_posxreal);
                                if moving < 0
                                    moving = moving+pi;
                                end
                            else
                                moving = 0;
                            end
                            maxangle = 1.5;
                            if cent2(1) > 1200
                                maxangle = 1.5;
                            end
                            delta_orientation = Delta_Angle(angle1, angle2);
                            
                            if delta_orientation < maxangle
                                
                                cvois = cvois+1;
                                lisvois(cvois) = l;
                                
                                area1 = boundarea1(single_cell_indices1(f1un(l)));
                                polyarea(bound1(:,2), bound1(:,1));
                            else
                                cvoisput = cvoisput+1;
                                lisvoisput(cvoisput) = l;
                            end
                        end
                    end
                end
                if cvois == 0
                    if cvoisput == 1
                        l = lisvoisput(1);
                        angle1 = slope1(single_cell_indices1(f1un(l)));
                        area1 = boundarea1(single_cell_indices1(f1un(l)));
                        long1 = longueur1(single_cell_indices1(f1un(l)));
                        if abs(long2-long1)/(long2+long1) < 0.2 && abs(area2-area1)/(area2+area1) < 0.1
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = 0;
                            f2untemp(l2,:) = 0;
                        end
                    elseif cvoisput == 2
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear distance
                        for k = 1:cvoisput
                            l = lisvoisput(k);
                            
                            bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                            centa{k} = center1(single_cell_indices1(f1un(l)),:);
                            anglea(k) = slope1(single_cell_indices1(f1un(l)));
                            longa(k) = longueur1(single_cell_indices1(f1un(l)));
                        end
                        if ((bounda(k)+bounda(k-1))/area2 > 0.9 && (bounda(k)+bounda(k-1))/area2 < 1.1)
                            meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                            if abs(meancent(1)- cent2(1)) < 10 + 1 && abs(meancent(2) - cent2(2)) < 10
                                for k = 1:cvoisput
                                    l = lisvoisput(k);
                                    pair(match_compt,:) = [l l2];
                                    match_compt = match_compt+1;
                                    f1untemp(l,:) = 0;
                                    f2untemp(l2,:) = 0;
                                    
                                end
                            end
                            
                        end
                        
                    end
                elseif cvois == 1
                    l = lisvois(1);
                    
                    
                    angle1 = slope1(single_cell_indices1(f1un(l)));
                    area1 = boundarea1(single_cell_indices1(f1un(l)));
                    long1 = longueur1(single_cell_indices1(f1un(l)));
                    delta_orientation = Delta_Angle(angle1, angle2);
                    if abs(long2-long1)/(long2+long1) < 0.2 && delta_orientation < 0.3 || abs(area2-area1)/(area2+area1) < 0.1 && delta_orientation < 0.3 || abs(area2-area1)/(area2+area1) < 0.25 && abs(long2-long1)/(long2+long1) < 0.3
                        pair(match_compt,:) = [l l2];
                        match_compt = match_compt+1;
                        f1untemp(l,:) = 0;
                        f2untemp(l2,:) = 0;
                    end
                    
                elseif cvois == 2
                    clear bounda;
                    clear centa;
                    clear anglea;
                    clear longa;
                    clear distance
                    for k = 1:cvois
                        l = lisvois(k);
                        
                        bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                        centa{k} = center1(single_cell_indices1(f1un(l)),:);
                        anglea(k) = slope1(single_cell_indices1(f1un(l)));
                        longa(k) = longueur1(single_cell_indices1(f1un(l)));
                        distance(k) = sqrt(sum(diff([centa{k};cent2]).^2));
                        disX(k) = centa{k}(1)-cent2(1);
                        disY(k) = centa{k}(2)-cent2(2);
                        
                    end
                    delta_orientation = Delta_Angle(anglea(k), anglea(k-1));
                    if (abs((sum(longa)-long2))/(sum(longa)+long2) < 0.1 && delta_orientation < 0.7) || ((sum(longa)-long2)/(sum(longa)+long2) < 0.25 && (bounda(k)+bounda(k-1))/area2 > 0.75 && (bounda(k)+bounda(k-1))/area2 < 1.25 && delta_orientation < 0.30 ) ||((sum(longa)-long2)/(sum(longa)+long2) < 0.2 && (bounda(k)+bounda(k-1))/area2 > 0.85 && (bounda(k)+bounda(k-1))/area2 < 1.15)
                        meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                        if abs(meancent(1)- cent2(1)) < max([mean_move 11]) && abs(meancent(2) - cent2(2)) < 10
                            for k = 1:cvois
                                l = lisvois(k);
                                pair(match_compt,:) = [l l2];
                                match_compt = match_compt+1;
                                f1untemp(l,:) = 0;
                                f2untemp(l2,:) = 0;
                                
                            end
                        end
                        
                    else
                        [valmindis idxmindis] = min(distance(:));
                        if abs(diff(distance)) < 1.5
                            [valmindis idxmindis] = min(abs(disY(:)));
                        end
                        delta_orientation = Delta_Angle(anglea(idxmindis), angle2);
                        if abs(longa(idxmindis) - long2)/(long2+longa(idxmindis)) < 0.1 && valmindis < 5 %&& delta_orientation < 0.4
                            l = lisvois(idxmindis);
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = 0;
                            f2untemp(l2,:) = 0;
                        end
                        
                    end
                elseif cvois > 2
                    cvois2 = 0;
                    clear lisvois2;
                    clear bounda;
                    clear centa;
                    clear anglea;
                    clear longa;
                    for k = 1:cvois
                        l = lisvois(k);
                        bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                        centa{k} = center1(single_cell_indices1(f1un(l)),:);
                        anglea(k) =slope1(single_cell_indices1(f1un(l)));
                        longa(k) = longueur1(single_cell_indices1(f1un(l)));
                        delta_posxreal(k) = centa{k}(1)-cent2(1);
                        delta_posyreal(k) = centa{k}(2)-cent2(2);
                        distance(k) = sqrt(sum(diff([centa{k};cent2]).^2));
                        if long2 < 40
                            if abs(delta_posxreal(k)) > 1 || abs(delta_posyreal(k)) > 1
                                moving(k) = atan(delta_posyreal(k)/delta_posxreal(k));
                                if moving(k) < 0
                                    moving(k) = moving(k)+pi;
                                end
                            else moving(k) = 0;
                            end
                        else
                            if abs(delta_posxreal(k)) > long2/30 || abs(delta_posyreal(k)) > long2/30;
                                moving(k) = atan(delta_posyreal(k)/delta_posxreal(k));
                                if moving(k) < 0
                                    moving(k) = moving(k)+pi;
                                end
                            else moving(k) = 0;
                            end
                        end
                        
                        if (Delta_Angle(moving(k), 0) < 0.2 || Delta_Angle(moving(k),angle2) < 0.7) && abs(delta_posyreal(k)) < max([sin(angle2)*long2 3.5])
                            cvois2 = cvois2+1;
                            lisvois2(cvois2) = l;
                            
                        end
                    end
                    
                    if cvois2 == 3 && long2 > 95
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear delta_posxreal;
                        clear delta_posyreal;
                        clear distance;
                        for k = 1:cvois2
                            l = lisvois2(k);
                            bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                            centa{k} = center1(single_cell_indices1(f1un(l)),:);
                            anglea(k) =slope1(single_cell_indices1(f1un(l)));
                            longa(k) = longueur1(single_cell_indices1(f1un(l)));
                            delta_posxreal(k) = centa{k}(1)-cent2(1);
                            delta_posyreal(k) = centa{k}(2)-cent2(2);
                            distance(k) = sqrt(sum(diff([centa{k};cent2]).^2));
                        end
                        [mintomaxlong, IX] = sort(longa);
                        for k = 2:3
                            l = lisvois2(IX(k));
                            pair(match_compt,:) = [l l2];
                            match_compt = match_compt+1;
                            f1untemp(l,:) = 0;
                            f2untemp(l2,:) = 0;
                        end
                        
                    elseif cvois2 == 3 && long2 < 95
                        cvois3 = 0;
                        clear lisvois3;
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear delta_posxreal;
                        clear delta_posyreal;
                        clear moving
                        clear distance
                        for k = 1:cvois2
                            l = lisvois2(k);
                            bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                            centa{k} = center1(single_cell_indices1(f1un(l)),:);
                            anglea(k) =slope1(single_cell_indices1(f1un(l)));
                            longa(k) = longueur1(single_cell_indices1(f1un(l)));
                            delta_posxreal(k) = centa{k}(1)-cent2(1);
                            delta_posyreal(k) = centa{k}(2)-cent2(2);
                            distance(k) = sqrt(sum(diff([centa{k};cent2]).^2));
                        end
                        bound1 = PerimPixel{Previous}{single_cell_indices1(f1un(lisvois2(1)))};
                        bound2 = PerimPixel{Previous}{single_cell_indices1(f1un(lisvois2(2)))};
                        bound3 = PerimPixel{Previous}{single_cell_indices1(f1un(lisvois2(3)))};
                        
                        [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound2(:,2)], [bound1(:,1); bound2(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(1) = first;
                            longue(1) = second;
                        else
                            width(1) = second;
                            longue(1) = first;
                        end
                        
                        
                        [rectx, recty, aire, perim] =  minboundrect([bound1(:,2); bound3(:,2)], [bound1(:,1); bound3(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(2) = first;
                            longueur(2) = second;
                        else
                            width(2) = second;
                            longue(2) = first;
                        end
                        
                        
                        [rectx, recty, aire, perim] =  minboundrect([bound3(:,2); bound2(:,2)], [bound3(:,1); bound2(:,1)]);
                        first = sqrt((rectx(2)-rectx(1))^2 + (recty(2)-recty(1))^2);
                        second = sqrt((rectx(3)-rectx(2))^2 + (recty(3)-recty(2))^2);
                        if first < second
                            width(3) = first;
                            longue(3) = second;
                        else
                            width(3) = second;
                            longue(3) = first;
                        end
                        [valmin minidx] = min(width);
                        k = minidx;
                        if k == 1
                            lisvois3 = [lisvois2(1) lisvois2(2)];
                        elseif k == 2
                            lisvois3 = [lisvois2(1) lisvois2(3)];
                        elseif k == 3
                            lisvois3 = [lisvois2(2) lisvois2(3)];
                        end
                        longueur = longue(k);
                        clear bounda;
                        clear centa;
                        clear anglea;
                        clear longa;
                        clear delta_posxreal;
                        clear delta_posyreal;
                        clear moving
                        clear distance
                        for k = 1:2
                            l = lisvois3(k);
                            bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                            centa{k} = center1(single_cell_indices1(f1un(l)),:);
                            anglea(k) =slope1(single_cell_indices1(f1un(l)));
                            longa(k) = longueur1(single_cell_indices1(f1un(l)));
                            delta_posxreal(k) = centa{k}(1)-cent2(1);
                            delta_posyreal(k) = centa{k}(2)-cent2(2);
                            distance(k) = sqrt(sum(diff([centa{k};cent2]).^2));
                        end
                        
                        delta_or = Delta_Angle(anglea(k), anglea(k-1));
                        if (delta_or < 0.3 && ((bounda(k)+bounda(k-1))/area2 > 0.75 && (bounda(k)+bounda(k-1))/area2 < 1.25) || abs((longueur-long2))/(longueur+long2) < 0.2 ) || (delta_or < 0.8 && (bounda(k)+bounda(k-1))/area2 > 0.85 && (bounda(k)+bounda(k-1))/area2 < 1.15)
                            meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                            if abs(meancent(1)- cent2(1)) < max([abs(cos(angle2))*long2/2, max_move]) && abs(meancent(2) - cent2(2)) <  abs(sin(angle2))*long2/2+1
                                for k = 1:2
                                    l = lisvois3(k);
                                    pair(match_compt,:) = [l l2];
                                    match_compt = match_compt+1;
                                    f1untemp(l,:) = 0;
                                    f2untemp(l2,:) = 0;
                                end
                            end
                        end
                    elseif cvois2 == 2
                        clear bounda
                        clear centa
                        clear anglea
                        clear longa
                        clear moving;
                        for k = 1:cvois2
                            l = lisvois2(k);
                            bounda(k) = boundarea1(single_cell_indices1(f1un(l)));
                            centa{k} = center1(single_cell_indices1(f1un(l)),:);
                            anglea(k) = slope1(single_cell_indices1(f1un(l)));
                            longa(k) = longueur1(single_cell_indices1(f1un(l)));
                        end
                        delta_orientation = Delta_Angle(anglea(k), anglea(k-1));
                        if( (bounda(k)+bounda(k-1))/area2 > 0.725 && (bounda(k)+bounda(k-1))/area2 < 1.225 && delta_orientation< 0.35 ) || ( (bounda(k)+bounda(k-1))/area2 > 0.85 && (bounda(k)+bounda(k-1))/area2 < 1.15 && delta_orientation< 0.8 )|| ( (bounda(k)+bounda(k-1))/area2 > 0.8 && (bounda(k)+bounda(k-1))/area2 < 1.2 && (longa(k)+longa(k-1))/long2 < 1.1 &&  (longa(k)+longa(k-1))/long2 > 0.9  )
                            
                            meancent = [mean([centa{k}(1) centa{k-1}(1)]) mean([centa{k}(2) centa{k-1}(2)])];
                            if abs(meancent(1)- cent2(1)) < abs(cos(angle2))*long2/2+mean_move && abs(meancent(2) - cent2(2)) <  max([abs(sin(angle2))*long2/2 3.5])
                                for k = 1:cvois2
                                    l = lisvois2(k);
                                    pair(match_compt,:) = [l l2];
                                    match_compt = match_compt+1;
                                    f1untemp(l,:) = 0;
                                    f2untemp(l2,:) = 0;
                                    
                                end
                            end
                        end
                    end
                elseif cvois2 == 1
                    delta_orientation = Delta_Angle(angle1, angle2);
                    if  abs(long2-long1)/(long2+long1) < 0.2 && delta_orientation < 0.2 && abs(area2-area1)/(area2+area1) < 0.1
                        pair(match_compt,:) = [l l2];
                        match_compt = match_compt+1;
                        f1untemp(l,:) = 0;
                        f2untemp(l2,:) = 0;
                    end
                end
                
                
            end
        end
    end
end


if exist('pair', 'var') == 1
    for k = 1:size(pair,1);
        if I(f2un(pair(k,2))) == 0;
            I(f2un(pair(k,2))) = f1un(pair(k,1));
        else
            I2(f2un(pair(k,2))) = f1un(pair(k,1));
        end
    end
end


assignin('base', 'I', I)

for i = 1:length(I)
    if I(i) ~= 0
        CurrentCenter(i,1:2) = [CX{Current}(i), CY{Current}(i)]
        PreviousCenter(I(i),1:2) = [CX{Previous}(I(i)), CY{Previous}(I(i))]

    end
end

assignin('base', 'I', I)

 %% Display Image
    
    %Display Image from combined objects
    if strcmp(Var.Figure.Display, 'on')
        FigNum = find((strcmp(Var.Figure.List, 'TrackObj')));
        figure(FigNum(CallNum))
        plot(CurrentCenter(:,1),CurrentCenter(:,2), '.b')
        hold on
        plot(PreviousCenter(:,1),PreviousCenter(:,2), '.r')
        axis([1 Var.Analysis.ImgSize(2) 1 Var.Analysis.ImgSize(1)])
        set(gca, 'YDir', 'reverse')
        
        for m = 1:size(CurrentCenter,1)
            text(CurrentCenter(m,1),CurrentCenter(m,2), num2str(m), 'FontSize', 16);
        end
        for m = 1:size(PreviousCenter,1)
            text(PreviousCenter(m,1),PreviousCenter(m,2), num2str(m),'FontAngle', 'italic', 'FontSize', 16); %
        end
        hold off
        
    end