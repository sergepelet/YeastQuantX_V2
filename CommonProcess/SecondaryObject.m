function Var = SecondaryObject(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end
Debug = 0;

%% Assigne Ref and input images
PrimObjName = Var.Analysis.SecPrimObj{CallNum};
ImageName = Var.Analysis.SecIntImage{CallNum};
SecObjName = Var.Analysis.SecObjOUT{CallNum};

NumPix = Var.Analysis.SecNumPix{CallNum};
Method = Var.Analysis.SecMethod{CallNum};

%% Get Images

PrimObj =  Var.Img.(PrimObjName);
IntImg=  Var.Img.(ImageName);
SecObj = zeros(size(PrimObj));
%Get Label for Primary object
PrimLabel = Var.Measurements.(PrimObjName).SegLabel;

if Debug
    figure(100)
    imagesc(PrimObj)
end

if strcmpi(Method, 'Border')
    %Make BW image
    BW = zeros(size(PrimObj));
    BW(PrimObj >0) = 1;
    %Get primeter
    Perim = bwperim(BW);
    %Grow Perimeter
    SE = strel('disk',NumPix);
    Perim = imdilate(Perim, SE);
    SecObj = Perim.*PrimObj;
elseif strcmpi(Method, 'Perimeter')
    %Grow inside and outside of cell!
    %Make BW image
    BW = zeros(size(PrimObj));
    BW(PrimObj >0) = 1;
    Dist = bwdist(BW);
    WatershedRegions = watershed(Dist);
    WatershedRegions(WatershedRegions >0) = 1;
    WatershedRegions = double (WatershedRegions);
    %NumPix inside
    PixIN = floor(NumPix);
    %NumPix outside given by first decimal
    PixOUT = round(rem(NumPix,1)*10);
    
    if PixIN > 0
        %Get perimeter
        Perim = bwperim(BW);
        %Grow Perimeter inside
        SE = strel('disk',PixIN);
        Perim = imdilate(Perim, SE);
        BorderIN = Perim.*PrimObj;
        SecObj = BorderIN;
    else
        SecObj = zeros(size(PrimObj));
    end
    if PixOUT >0
        %Loop through all objects
        for O = 1:max(PrimObj(:))
            %Create empty images which will contain only one object
            CenterObj = zeros(size(PrimObj));
            %Place Small Object in CentreObj
            CenterObj(PrimObj == O) = 1;
            %Greate extended perimeter object outside of CenterObj
            SE = strel('disk',PixOUT);
            BorderOUT = imdilate(CenterObj, SE)-CenterObj;
            %Find Watershed region belonging to object
            
            SingleWR = imextendedmax(WatershedRegions+CenterObj,1);
            
            %Remove Pix that do not belong to WS region
            BorderOUT = BorderOUT.*SingleWR;
            %Add Borderout Pix to secondary object image
            SecObj(BorderOUT ==1) = O;
        end
    end
elseif strcmpi(Method, 'Center')
    %Make BW image
    BW = zeros(size(PrimObj));
    BW(PrimObj >0) = 1;
    %Shrink objects to points
    SecObj = bwmorph(BW,'shrink',Inf);
    %Grow objects to desired size
    SE = strel('disk',NumPix);
    SecObj = imdilate(SecObj, SE);
    %Remove pixels outside of primary object
    SecObj = SecObj.*PrimObj;
    
    
elseif strcmpi(Method, 'TMD')
    %Make BW image
    BW = zeros(size(PrimObj));
    BW(PrimObj >0) = 1;
    %Get perimeter
    Perim = bwperim(BW);
    %Grow Perimeter
    SE = strel('disk',NumPix);
    Perim = imdilate(Perim, SE);
    MembraneObj = Perim.*PrimObj;
    [MembraneObj , NumObj] = bwlabel(MembraneObj);
    
    ObjProps = regionprops(MembraneObj,IntImg,'PixelIdxList','PixelValues','PixelList','Centroid');
    
    
    
    HiPixMembrane = zeros(size(PrimObj));
    
    
    for O = 1:NumObj
        %find Slices coordinates
        Qind = SplitSlices(ObjProps(O).PixelList(:,1),ObjProps(O).PixelList(:,2), ObjProps(O).Centroid(1), ObjProps(O).Centroid(2));
        
        NumMemPix = length(ObjProps(O).PixelIdxList);
        %Find hi pixels in each slice
        %NumSelect = round(NumPix/16);
        NumSelect = round(NumMemPix/NumPix*1.5/16);
        for S = 1:16
            [SortInt, SortInd] = sort(ObjProps(O).PixelValues(Qind{S}));
            if length (SortInd) > NumSelect
                PixIn = SortInd(end-NumSelect:end);
            else
                PixIn = SortInd(:);
            end
            HiPixMembrane(ObjProps(O).PixelIdxList(Qind{S}(PixIn))) = 1;
        end
        
    end
    
    SE = strel('disk',2);
    SecObj = imclose(HiPixMembrane, SE);
    SecObj = SecObj.*PrimObj;
    %
    %         figure(100)
    %         subplot(2,2,1); imagesc(IntImg);title('IntImg')
    %         subplot(2,2,2); imagesc(MembraneObj); title('MembraneObj');
    %         subplot(2,2,3); imagesc(HiPixMembrane); title('HiPixMembrane');
    %         AllObj = zeros(size(PrimObj));
    %         AllObj(PrimObj >0) = 1;
    %        % AllObj(MembraneObj >0) = AllObj(MembraneObj >0) + 1;
    %        % AllObj(HiPixMembrane >0) = AllObj(HiPixMembrane >0) + 2;
    %         AllObj(SecObj >0) = AllObj(SecObj >0) + 2;
    %         subplot(2,2,4); imagesc(AllObj); title('AllObj');
elseif strcmpi(Method, 'MaxPixels') ||strcmpi(Method, 'MinPixels')
    
    
    
    %% Get Object and Image Properties
    ObjProps = regionprops(PrimObj,IntImg,'PixelIdxList','PixelValues','Centroid');
    for O = 1:max(PrimObj(:))
        [SortInt, SortInd] = sort(ObjProps(O).PixelValues);
        if strcmpi(Method, 'MaxPixels')
            if length (SortInd) > NumPix
                PixIn = SortInd(end-NumPix:end);
            else
                PixIn = SortInd(:);
            end
        else strcmpi(Method, 'MinPixels')
            if length (SortInd) < NumPix
                PixIn = SortInd (:);
            else
                PixIn = SortInd(1:NumPix);
            end
        end
        SecObj(ObjProps(O).PixelIdxList(PixIn)) = O;
        
        
        
        Var.Measurements.(SecObjName).CenterX(PrimLabel(O),1) = ObjProps(O).Centroid(1);
        Var.Measurements.(SecObjName).CenterY(PrimLabel(O),1) = ObjProps(O).Centroid(2);
        Var.Measurements.(SecObjName).PixelList{1,PrimLabel(O)} = PixIn;
        
        
    end
elseif strcmpi(Method, 'MaxPixelsFraction') ||strcmpi(Method, 'MinPixelsFraction')
    
    
    
    %% Get Object and Image Properties
    ObjProps = regionprops(PrimObj,IntImg,'PixelIdxList','PixelValues','Centroid');
    
    
    for O = 1:max(PrimObj(:))
      % O=O
       
           
        [SortInt, SortInd] = sort(ObjProps(O).PixelValues);
       SizeObj = ceil(length(SortInt));
        NumberPixels = ceil(length(SortInt)*NumPix/100);
        if strcmpi(Method, 'MaxPixelsFraction')
            PixIn = SortInd(end-NumberPixels:end);
        else strcmpi(Method, 'MinPixelsFraction')
            PixIn = SortInd(1:NumberPixels);
        end
        SecObj(ObjProps(O).PixelIdxList(PixIn)) = O;
        
        
        
        Var.Measurements.(SecObjName).CenterX(PrimLabel(O),1) = ObjProps(O).Centroid(1);
        Var.Measurements.(SecObjName).CenterY(PrimLabel(O),1) = ObjProps(O).Centroid(2);
        Var.Measurements.(SecObjName).PixelList{1,PrimLabel(O)} = PixIn;
        
        
    end
elseif(strcmpi(Method, 'AboveThreshold'))
    %% Get Object and Image Properties
    ObjProps = regionprops(PrimObj,IntImg,'PixelIdxList','PixelValues','Centroid');
    for O = 1:max(PrimObj(:))
        %Find pixels above the threshold
        PixAbove = find(ObjProps(O).PixelValues >NumPix);
        SecObj(ObjProps(O).PixelIdxList(PixAbove)) = O;
        
        Var.Measurements.(SecObjName).CenterX(PrimLabel(O),1) = ObjProps(O).Centroid(1);
        Var.Measurements.(SecObjName).CenterY(PrimLabel(O),1) = ObjProps(O).Centroid(2);
        Var.Measurements.(SecObjName).PixelList{1,PrimLabel(O)} = PixAbove;
    end
    
elseif(strcmpi(Method, 'BelowThreshold'))
    %% Get Object and Image Properties
    ObjProps = regionprops(PrimObj,IntImg,'PixelIdxList','PixelValues','Centroid');
    for O = 1:max(PrimObj(:))
        %Find pixels above the threshold
        PixBelow = find(ObjProps(O).PixelValues <NumPix);
        SecObj(ObjProps(O).PixelIdxList(PixBelow)) = O;
        
        Var.Measurements.(SecObjName).CenterX(PrimLabel(O),1) = ObjProps(O).Centroid(1);
        Var.Measurements.(SecObjName).CenterY(PrimLabel(O),1) = ObjProps(O).Centroid(2);
        Var.Measurements.(SecObjName).PixelList{1,PrimLabel(O)} = PixBelow;
    end
    %%
elseif(strcmpi(Method, 'Circle'))
    %% Get Object and Image Properties
    ObjProps = regionprops(PrimObj,IntImg,'PixelIdxList','Centroid', 'EquivDiameter', 'MinorAxisLength');
    SecObj = zeros(size(PrimObj));
    %                 figure(100)
    %     imagesc(PrimObj)
    %     hold on
    
    
    %     OUTImg = zeros(size(PrimObj));
    %     OUTImg(PrimObj >=1) =1;
    
    % generate empty image
    BW = zeros(size(PrimObj));
    Angle = [0 : (2 * pi / 100) : (2 * pi)];
    for O = 1:max(PrimObj(:))
        %Set object in image to one
        BW(:) = 0;
        BW(ObjProps(O).PixelIdxList) = 1;
        
        %get coordinates from pixels in object
        [YGr,XGr] = find(BW > 0);
        
        %Find circle in image
        [CellCenter,CellRad] = imfindcircles(BW,[floor(ObjProps(O).EquivDiameter/4),ceil(ObjProps(O).EquivDiameter/1.8)],'Sensitivity',0.95);
        
        
        %Sort Radius to start with larger circle
        [CellRad, RadOrder] = sort(CellRad, 'descend');
        CellCenter = CellCenter(RadOrder,:);
        
        % if no circle found generate circle from centroid and MinorAxis
        if isempty(CellRad)
            CellRad(1) = ObjProps(O).MinorAxisLength/2;
            CellCenter = ObjProps(O).Centroid;
            %             viscircles(CellCenter, CellRad,'EdgeColor','r');
            %         pause(0.1)
       % else
            %             viscircles(CellCenter(1,:), CellRad(1),'EdgeColor','k');
            %             pause(0.1)
        end
        
        % Use the largest circle
        % find pixels within circle given by radius and center
        
        for C = 1:1 %length(CellRad)
            XPerim = CellRad(C) * cos(Angle) + CellCenter(C,1);
            YPerim = CellRad(C) * sin(Angle) + CellCenter(C,2);
            InCell = inpolygon(XGr,YGr,XPerim,YPerim);
            %Find indices within circle
            InInd = sub2ind(size(BW),YGr(InCell),XGr(InCell));
            %             OUTImg(InInd) = C+3;
            %
            %             figure(100)
            %             imagesc(OUTImg)
            %pause(1)
        end
        
        
        
        %Set secondary image to pixels in the cell
        SecObj(InInd) = O;
        %Set measurements data for center and pixel list
        Var.Measurements.(SecObjName).CenterX(PrimLabel(O),1) = CellCenter(1,1);
        Var.Measurements.(SecObjName).CenterY(PrimLabel(O),1) = CellCenter(1,2);
        Var.Measurements.(SecObjName).PixelList{1,PrimLabel(O)} = InInd;
    end
end


Var.Measurements.(SecObjName).SegLabel = PrimLabel;




%% Display %%%
if strcmp(Var.Figure.Display, 'on')
    FigNum = find(strcmp(Var.Figure.List, 'SecPrimObj'));
    figure(FigNum(CallNum))
    subplot(2,2,1); imagesc(PrimObj);title('Primary Object Image')
    subplot(2,2,2); imagesc(IntImg); title('Intenisty Image');
    subplot(2,2,3); imagesc(SecObj); title('Secondary Object');
    Obj = zeros(size(SecObj));
    Obj(SecObj>0) = 1;
    BothObjRGB(:,:,1) = Obj;
    Obj = zeros(size(PrimObj));
    Obj(PrimObj>0) = 0.5;
    BothObjRGB(:,:,3) = Obj;
    BothObjRGB(:,:,2) = (IntImg-min(IntImg(:)))./(max(IntImg(:))-min(IntImg(:)));
    subplot(2,2,4); image(BothObjRGB); title('Overlay');
    
    
end


%%% save %%%
Var.Img.(SecObjName) = SecObj;

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%




function Qind = SplitSlices(Xpix,Ypix, X0, Y0)

%Split object in 16 slices
%Qind =  the indices of the split slices

% Xpix = ObjProps(O).PixelList(:,1);
% Ypix = ObjProps(O).PixelList(:,2);
% X0 = ObjProps(O).Centroid(1);
% Y0 = ObjProps(O).Centroid(2);
for Q = 1:16
    switch Q
        case 1
            Qind{Q} = intersect( find((Xpix > X0)), find(Ypix > 2*Xpix- 2*X0 +Y0));
            
        case 2
            Qind{Q} = intersect(find(Ypix < 2*Xpix- 2*X0 +Y0),  find(Ypix > Xpix- X0 +Y0));
            
        case 3
            Qind{Q} = intersect(find(Ypix < Xpix- X0 +Y0),  find(Ypix > 0.5*Xpix- 0.5*X0 +Y0));
            
        case 4
            Qind{Q} = intersect( find(Ypix < 0.5*Xpix- 0.5*X0 +Y0), find(Ypix >Y0));
            
        case 5
            Qind{Q} = intersect(find(Ypix > -0.5*Xpix + 0.5*X0 + Y0), find(Ypix < Y0));
            
        case 6
            Qind{Q} = intersect(find(Ypix > -1*Xpix + X0 + Y0),find(Ypix < -0.5*Xpix + 0.5*X0 + Y0));
            
        case 7
            Qind{Q} = intersect(find(Ypix > -2*Xpix + X0*2 + Y0),find(Ypix < -1*Xpix + X0 + Y0));
            
        case 8
            Qind{Q} = intersect(find((Xpix > X0)) , find(Ypix < -2*Xpix + X0*2 + Y0));
            
        case 9
            Qind{Q} = intersect( find((Xpix < X0)), find(Ypix < 2*Xpix- 2*X0 +Y0));
            
        case 10
            Qind{Q} = intersect(find(Ypix > 2*Xpix- 2*X0 +Y0),  find(Ypix < Xpix- X0 +Y0));
            
        case 11
            Qind{Q} = intersect(find(Ypix > Xpix- X0 +Y0),  find(Ypix < 0.5*Xpix- 0.5*X0 +Y0));
            
        case 12
            Qind{Q} = intersect( find(Ypix > 0.5*Xpix- 0.5*X0 +Y0), find(Ypix <Y0));
            
        case 13
            Qind{Q} = intersect(find(Ypix < -0.5*Xpix + 0.5*X0 + Y0), find(Ypix > Y0));
            
        case 14
            Qind{Q} = intersect(find(Ypix < -1*Xpix + X0 + Y0),find(Ypix > -0.5*Xpix + 0.5*X0 + Y0));
            
        case 15
            Qind{Q} = intersect(find(Ypix < -2*Xpix + X0*2 + Y0),find(Ypix > -1*Xpix + X0 + Y0));
            
        case 16
            Qind{Q} = intersect(find((Xpix < X0)) , find(Ypix > -2*Xpix + X0*2 + Y0));
            
    end
end
