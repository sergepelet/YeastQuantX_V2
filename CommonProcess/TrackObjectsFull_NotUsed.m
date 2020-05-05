function Var = TrackObjectsFull(Var, CallNum)
tic
Debug = 0;      %Set Debug to 0 to prevent image display

if nargin == 1
    CallNum = 1;
end
%Get Track Object and Tracking distance from Var
TrackObj = 'Nucl' ; Var.Analysis.TrackObj{CallNum};
MaxDistance = Var.Analysis.TrackMaxDistance{CallNum};

%Initalize images
%TrackImg = zeros(Var.Analysis.ImgSize);
TrackCenterSum = zeros(Var.Analysis.ImgSize);
TrackSum = zeros(Var.Analysis.ImgSize);

%Loop through all frames
for F = 1:Var.Analysis.NumFrame
    TrackImg = zeros(Var.Analysis.ImgSize);
    %Read Frame.mat file
    ImgFolder = fileparts(Var.Analysis.OutPath);
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(F,'%05d'),'.mat'])
    load(FileName)
    
    %Get X and Y centers from measurement
    CX{F} = round(Var.Measurements.(TrackObj).CenterX);
    CY{F} = round(Var.Measurements.(TrackObj).CenterY);
    
    %Place centers of object in image
    for Obj = 1:length(CX{F})
        TrackImg(CY{F}(Obj),CX{F}(Obj)) = 1;
    end
    
    %Dilate centers and add them in a single image
    SE = strel('disk', MaxDistance/2);
    TrackArea = imdilate(TrackImg, SE);
    TrackSum = TrackSum + TrackArea;
    
    TrackCenterSum = TrackCenterSum + TrackImg;
    
end
figure(100)
imagesc(TrackImg)

figure(102)
imagesc(TrackSum)
pause(0.5)
%Create image which contains only the objects tracked from begining to end
TrackTop = zeros(Var.Analysis.ImgSize);
TrackTop(TrackSum == Var.Analysis.NumFrame) = 1;

DisplayTrackTop = bwlabel(TrackTop);
DisplayTrackTop = DisplayTrackTop + 10.*TrackCenterSum;
figure(103)
imagesc(DisplayTrackTop)
pause(0.5)

%Get centroids from these objects
CC = bwconncomp(TrackTop);
Centro = regionprops(CC,'centroid');

assignin('base','Centro',Centro);

%Loop trhough all centroids
for C = 1:length(Centro)
    %For each frame find corresponding Nucleus
    for F = 1:Var.Analysis.NumFrame
        if C== 1
            CenterList{F} = [];
        end
        GoodX = find(CX{F}>=Centro(C).Centroid(1)-MaxDistance & CX{F}<=Centro(C).Centroid(1)+MaxDistance);
        GoodY = find(CY{F}>=Centro(C).Centroid(2)-MaxDistance & CY{F}<=Centro(C).Centroid(2)+MaxDistance);
        GoodXY = intersect(GoodX,GoodY);
        if length(GoodXY) == 1
            GoodCenter{C,F} = GoodXY;
            CenterList{F} = [CenterList{F}; GoodXY ];
        elseif length(GoodXY) > 1
            XYDist = [];
            for X = 1:length(GoodXY)
                XYDist(X) = sqrt((Centro(C).Centroid(1)-CX{F}(GoodXY(X)))^2 + (Centro(C).Centroid(2)-CY{F}(GoodXY(X)))^2);
            end
            [XYDist, SortInd] = sort(XYDist, 'ascend');
            GoodXY = GoodXY(SortInd);
            
            BelowMax = find(XYDist <= MaxDistance,1, 'last');
            if ~isempty(BelowMax)
                
                GoodCenter{C,F} = GoodXY(1:BelowMax);
                CenterDist{C,F} = XYDist(1:BelowMax);
                CenterList{F} = [CenterList{F}; GoodXY(1:BelowMax) ];
            end
            
        end
    end
end
assignin('base','GoodCenter',GoodCenter);
assignin('base','CenterDist',CenterDist);


NumCenter = cellfun(@length,GoodCenter);
for C = 1:length(Centro)
    if max(NumCenter(C,:)) == 1 && min(NumCenter(C,:)) == 1
        for F = 1:Var.Analysis.NumFrame
            NumObj = GoodCenter{C,F};
            IndObj = find(CenterList{F}(:) == NumObj);
            if length(IndObj)== 1
                TrackMat(C,F) = NumObj;
            elseif IndObj > 1
                F = F
                NumObj = NumObj
                IndObj = IndObj
                
                find
                
            end
        end
    end
end

assignin('base','TrackMat',TrackMat);

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;




