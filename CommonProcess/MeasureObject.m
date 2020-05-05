function Var = MeasureObject(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end

ImgName = Var.Analysis.MeasureIntensity{floor(CallNum)};
ObjName = Var.Analysis.MeasureObject{round(mod(CallNum,1)*100)};

CurrentFrame = Var.Analysis.CurrentFrame;

ObjectImage = Var.Img.(ObjName);
IntensityImage = Var.Img.(ImgName);

% %Bit of code to create an artificial image to load with a Data.mat file
% % that does not contain the images
% IntensityImage = zeros(Var.Analysis.ImgSize);
% ObjectImage = zeros(Var.Analysis.ImgSize);
% PixL = Var.Measurements.(ObjName).PixelList;
% PixCell = Var.Measurements.Cell.PixelList;
% T = 34;
% iter = 0;
% for O = 1:size(PixL, 2)
%     if ~isempty(PixL)
%         iter = iter+1;
%         ObjectImage(PixL{T,O}) = 1;
%         IntensityImage(PixCell{T,O}) = iter + rand(size(PixCell{T,O}));
%     end
% end
% ObjectImage = bwlabel(ObjectImage);

if isfield(Var.Analysis, 'NumHiPix')
    NumHiPix = Var.Analysis.NumHiPix;
else
    NumHiPix = 20;
end
% figure(100)
% imagesc(ObjectImage)

ObjectsProps = regionprops(ObjectImage,IntensityImage, 'Area', 'PixelIdxList','PixelList','PixelValues' , 'Centroid', ...
    'EquivDiameter','Orientation','MajorAxisLength','MinorAxisLength', 'Eccentricity','ConvexArea','MeanIntensity','WeightedCentroid');
NumObjects = size(ObjectsProps,1);


%If polarity has to be measured
if isfield(Var.Analysis, 'PolAnalysis') && strcmp(Var.Analysis.PolAnalysis, 'Yes')
    %get Database info
    DegResolution = Var.Analysis.PolDegree;
    Orientation = Var.Analysis.PolOrient;
    
    %generate list of angles
    AngleList = [-180:DegResolution:180];
    Var.Analysis.AngleAxis = AngleList(1:end-1)+ceil(DegResolution/2);
    %Generage Blank Image
    PolImage = zeros(size(ObjectImage));
    DistImage = zeros(size(ObjectImage));
    %Get pixelList for perimeters of all objects
    Perim = ObjectImage;
    Perim(ObjectImage>0) = 1;
    Perim = bwperim(Perim);
    PerimPix = find(Perim >0);
end


if isfield(Var, 'Measurements')  && isfield(Var.Measurements, ObjName)&& isfield(Var.Measurements.(ObjName), 'SegLabel')
    Label = Var.Measurements.(ObjName).SegLabel;
    for O = 1:NumObjects
        %Object Present in Frame
        %Size Features
        Var.Measurements.(ObjName).PixelList{1,O} = ObjectsProps(O).PixelIdxList;
        Var.Measurements.(ObjName).(ImgName).Area(O) = ObjectsProps(O).Area;
       % A = ObjectsProps(O).Area
        Var.Measurements.(ObjName).(ImgName).Diameter(O) = ObjectsProps(O).EquivDiameter;
        Var.Measurements.(ObjName).(ImgName).CenterX(O) = ObjectsProps(O).Centroid(1);
        Var.Measurements.(ObjName).(ImgName).CenterY(O) = ObjectsProps(O).Centroid(2);
        
        Var.Measurements.(ObjName).(ImgName).CheckLabel(O) = Label(O);
        %Get label of parent object if it exists
        if isfield(Var.Measurements.(ObjName), 'ParentLabel')
            Var.Measurements.(ObjName).(ImgName).ParentLabel(O) = Var.Measurements.(ObjName).ParentLabel(O);
        end
        
        %Intensity features
        Var.Measurements.(ObjName).(ImgName).TotalIntensity(O) = sum(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).AverageIntensity(O) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MedianIntensity(O) = median(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).STDIntensity(O) = std(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MaxIntensity(O) = max(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MinIntensity(O) = min(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).STDIntensityNorm(O) = std(IntensityImage(ObjectsProps(O).PixelIdxList))./mean(IntensityImage(ObjectsProps(O).PixelIdxList));
        %HiPix calculation
        if ObjectsProps(O).Area <= NumHiPix
            Var.Measurements.(ObjName).(ImgName).LoPix(O) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            Var.Measurements.(ObjName).(ImgName).HiPix(O) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            Var.Measurements.(ObjName).(ImgName).HiPix05(O) =  mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            Var.Measurements.(ObjName).(ImgName).HiPix10(O) =  mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            
            Var.Measurements.(ObjName).(ImgName).HiPixArea(O) = ObjectsProps(O).Area;
            Var.Measurements.(ObjName).(ImgName).HiPixList{1,O} = ObjectsProps(O).PixelIdxList;
            Var.Measurements.(ObjName).(ImgName).LoPixList{1,O} = ObjectsProps(O).PixelIdxList;
             Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(O) = NaN;
             Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(O) = NaN;
             Var.Measurements.(ObjName).(ImgName).ConnectedHiPixList{1,O} = [];
        else
            [SortedInt, SortedIndex] = sort(IntensityImage(ObjectsProps(O).PixelIdxList));
         
            Var.Measurements.(ObjName).(ImgName).LoPix(O) = mean(SortedInt(1:NumHiPix));
            Var.Measurements.(ObjName).(ImgName).HiPix(O) = mean(SortedInt(end-NumHiPix+1:end));
            
            Var.Measurements.(ObjName).(ImgName).HiPixList{1,O} = SortedIndex(end-NumHiPix+1:end);
            Var.Measurements.(ObjName).(ImgName).LoPixList{1,O} = SortedIndex(1:NumHiPix);
            
            %Various HiPix
            Var.Measurements.(ObjName).(ImgName).HiPix05(O) = mean(SortedInt(end-5+1:end));
            Var.Measurements.(ObjName).(ImgName).HiPix10(O) = mean(SortedInt(end-10+1:end));
            %HiPix Area
            HiPixImg = zeros(size(ObjectImage));
            HiPixImg(ObjectsProps(O).PixelIdxList(SortedIndex(end-NumHiPix+1:end))) = 1;
            HiPixStat = regionprops(HiPixImg, 'ConvexArea');
            Var.Measurements.(ObjName).(ImgName).HiPixArea(O) = HiPixStat.ConvexArea;
%             figure(103)
%             imagesc(ObjectImage +5.*HiPixImg)
            %Connected HiPix
            %Open image to remove isolated pixels
            SE = strel('disk',1);
            HiPixImg = imopen(HiPixImg,SE);
            %Label opened image
            [LabelHP, NumArea] = bwlabel(HiPixImg, 8);
%             figure(104)
%             imagesc(LabelHP)
%             pause(0.5)
            
            if NumArea == 1 %If a single area is generate just measure it
                ConnectedHiPixStat = regionprops(LabelHP, IntensityImage, 'Area', 'MeanIntensity', 'PixelIdxList');  %,
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(O) = ConnectedHiPixStat.Area;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(O) = ConnectedHiPixStat.MeanIntensity;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixList{1,O} = ConnectedHiPixStat.PixelIdxList;
            elseif NumArea == 0  %If no pixels are left after image opening set meas to NaN
                 Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(O) = NaN;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(O) = NaN;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixList{1,O} = [];
            else   %If multiple objects are generated find the biggest one and measure its intensity
                CHPArea = 0;
                CHPInt = 0;
                CHPList = [];
                for HPO = 1:NumArea
                    IndHPO = find(LabelHP == HPO);
                    if length(IndHPO) > CHPArea
                        CHPArea = length(IndHPO);
                        CHPInt = mean(IntensityImage(IndHPO));
                        CHPList = IndHPO;
                    end
                end
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(O) = CHPArea;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(O) = CHPInt;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixList{1,O} = CHPList;
            end
            
            
        end
        %Polarity Measurements using object shape
        Var.Measurements.(ObjName).(ImgName).Orientation(O) = 90-ObjectsProps(O).Orientation;
        Var.Measurements.(ObjName).(ImgName).MajorAxisLength(O) = ObjectsProps(O).MajorAxisLength;
        Var.Measurements.(ObjName).(ImgName).MinorAxisLength(O) = ObjectsProps(O).MinorAxisLength;
        Var.Measurements.(ObjName).(ImgName).Eccentricity(O) = ObjectsProps(O).Eccentricity;
        Var.Measurements.(ObjName).(ImgName).ConvexRatio(O) = ObjectsProps(O).Area/ObjectsProps(O).ConvexArea;
        
        %Polarity Measurements using weighted centroid
        Var.Measurements.(ObjName).(ImgName).WeightCenterX(O) = ObjectsProps(O).WeightedCentroid(1);
        Var.Measurements.(ObjName).(ImgName).WeightCenterY(O) = ObjectsProps(O).WeightedCentroid(2);
        WeightAxis = (ObjectsProps(O).WeightedCentroid(1)-ObjectsProps(O).Centroid(1))^2 + (ObjectsProps(O).WeightedCentroid(2)-ObjectsProps(O).Centroid(2))^2;
        Var.Measurements.(ObjName).(ImgName).WeightAxis(O) = WeightAxis;
        WeightAngle = atand(ObjectsProps(O).WeightedCentroid(1)-ObjectsProps(O).Centroid(1))/(ObjectsProps(O).WeightedCentroid(2)-ObjectsProps(O).Centroid(2));
        Var.Measurements.(ObjName).(ImgName).WeightAngle(O) = WeightAngle;
        
        %Polarity based on Polar coordinates
        if isfield(Var.Analysis, 'PolAnalysis') && strcmp(Var.Analysis.PolAnalysis, 'Yes')

            
            %Get pixels relative to center of mass
            Xpix = ObjectsProps(O).PixelList(:,1) - ObjectsProps(O).Centroid(1);
            Ypix = ObjectsProps(O).PixelList(:,2) - ObjectsProps(O).Centroid(2);
            %Convert to polar coordinates
            [Theta, Rho] = cart2pol(Xpix,Ypix);
            %Convert to degrees
            %And Rotate Angle to match desired orientation
             Angle = radtodeg(Theta);
            switch Orientation
                %         case '3h'
                %             Angle = radtodeg(Theta);
                case '6h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= -90) = AngleOrig(AngleOrig >= -90) -90;
                    Angle(AngleOrig < -90) = AngleOrig(AngleOrig < -90)+270;
                case '9h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= 0) = AngleOrig(AngleOrig >= 0) -180;
                    Angle(AngleOrig < 0) = 180 + AngleOrig(AngleOrig < 0);
                case '12h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= 90) = AngleOrig(AngleOrig >= 90) -270;
                    Angle(AngleOrig < 90) = AngleOrig(AngleOrig < 90)+90;
            end
            
            PerimInt = zeros(1,length(AngleList)-1);
            Distance2Center = zeros(1,length(AngleList)-1);
            PerimSize = zeros(1,length(AngleList)-1);
            for a = 1:length(AngleList)-1
                %find pixels in right angle range
                InRange = find(Angle > AngleList(a) & Angle <=AngleList(a+1));
                PerimInt(a) = mean(IntensityImage(ObjectsProps(O).PixelIdxList(InRange)));
                %get mean distance from center of all pixels in angle range 
                Distance2Center(a) = mean(Rho(InRange));
                
                %Calculate the intersect between perimeter and slice pixels
                PerimSize(a) = length(intersect(PerimPix, ObjectsProps(O).PixelIdxList(InRange)));
                
                %add pixels to image
                PolImage(ObjectsProps(O).PixelIdxList(InRange)) = a/length(AngleList);
                DistImage(ObjectsProps(O).PixelIdxList(InRange)) = Distance2Center(a);
            end
            %%% Display %%%
            if 0
                figure(101)
                subplot(2,2,1); imagesc(IntensityImage);title('IntensityImage')
                subplot(2,2,2); imagesc(ObjectImage); title('ObjectImage');
                subplot(2,2,3); imagesc(PolImage); title('Polarity image');
                subplot(2,2,4); imagesc(DistImage); title('Distance image');
                %pause(0.5)
            end
            
            %Add Measurement
            Var.Measurements.(ObjName).(ImgName).Polarity(O, 1:length(AngleList)-1 ) = PerimInt;

            [Var.Measurements.(ObjName).(ImgName).PolarityMaxInt(O), MaxInd] = max(PerimInt);
            Var.Measurements.(ObjName).(ImgName).PolarityMaxAngle(O) = Var.Analysis.AngleAxis(MaxInd);
            Var.Measurements.(ObjName).(ImgName).PolarityStdInt(O) = std(PerimInt);
            
            Var.Measurements.(ObjName).(ImgName).PerimSize(O, 1:length(AngleList)-1 ) = PerimSize;
            
            Var.Measurements.(ObjName).(ImgName).Distance2Center(O, 1:length(AngleList)-1 ) = Distance2Center;
            [Var.Measurements.(ObjName).(ImgName).D2C_Max(O), MaxInd] = max(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_MaxAngle(O) = Var.Analysis.AngleAxis(MaxInd);
            Var.Measurements.(ObjName).(ImgName).D2C_Mean(O) = mean(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_Std(O) = std(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_StdNorm(O) = std(Distance2Center)/mean(Distance2Center);
            
            
            %% Use Circle object as center for polarity measurement
            CircleMethod = find(strcmp(Var.Analysis.SecMethod, 'Circle'));
            
            if ~isempty(CircleMethod)
            CircleObjName = Var.Analysis.SecObjOUT{CircleMethod};
                      
            %Get pixels relative to center of mass
            Xpix = ObjectsProps(O).PixelList(:,1) - Var.Measurements.(CircleObjName).CenterX(O);
            Ypix = ObjectsProps(O).PixelList(:,2) - Var.Measurements.(CircleObjName).CenterY(O);
            
            
            %Convert to polar coordinates
            [Theta, Rho] = cart2pol(Xpix,Ypix);
            %Convert to degrees
            %And Rotate Angle to match desired orientation
             Angle = radtodeg(Theta);
            switch Orientation
                %         case '3h'
                %             Angle = radtodeg(Theta);
                case '6h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= -90) = AngleOrig(AngleOrig >= -90) -90;
                    Angle(AngleOrig < -90) = AngleOrig(AngleOrig < -90)+270;
                case '9h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= 0) = AngleOrig(AngleOrig >= 0) -180;
                    Angle(AngleOrig < 0) = 180 + AngleOrig(AngleOrig < 0);
                case '12h'
                    AngleOrig = radtodeg(Theta);
                    Angle(AngleOrig >= 90) = AngleOrig(AngleOrig >= 90) -270;
                    Angle(AngleOrig < 90) = AngleOrig(AngleOrig < 90)+90;
            end
            
            CC_PerimInt = zeros(1,length(AngleList)-1);
            CC_Distance2Center = zeros(1,length(AngleList)-1);
            for a = 1:length(AngleList)-1
                %find pixels in right angle range
                InRange = find(Angle > AngleList(a) & Angle <=AngleList(a+1));
                CC_PerimInt(a) = mean(IntensityImage(ObjectsProps(O).PixelIdxList(InRange)));
                %get mean distance from center of all pixels in angle range 
                CC_Distance2Center(a) = mean(Rho(InRange));
                
                %add pixels to image
                PolImage(ObjectsProps(O).PixelIdxList(InRange)) = a/length(AngleList);
                DistImage(ObjectsProps(O).PixelIdxList(InRange)) = CC_Distance2Center(a);
            end
            %%% Display %%%
            if 0
                figure(101)
                subplot(2,2,1); imagesc(IntensityImage);title('IntensityImage')
                subplot(2,2,2); imagesc(ObjectImage); title('ObjectImage');
                subplot(2,2,3); imagesc(PolImage); title('Polarity image');
                subplot(2,2,4); imagesc(DistImage); title('Distance image');
                %pause(0.5)
            end
            
            %Add Measurement
            Var.Measurements.(ObjName).(ImgName).Polarity2CC(O, 1:length(AngleList)-1 ) = CC_PerimInt;
            
            Var.Measurements.(ObjName).(ImgName).Distance2CC(O, 1:length(AngleList)-1 ) = CC_Distance2Center;
             
            end
            
            
        end
        
        %If Crowding has to be measured
        if isfield(Var.Analysis, 'CrowdingAnalysis') && strcmp(Var.Analysis.CrowdingAnalysis, 'Yes')
           
            CX = Var.Measurements.(ObjName).(ImgName).CenterX(O); 
            CY = Var.Measurements.(ObjName).(ImgName).CenterY(O);
            Distance = zeros(1,NumObjects);
            %Calculate distance between present object and all objects
            for P = 1:NumObjects
                Distance(P) = sqrt((CX- ObjectsProps(P).Centroid(1))^2+(CY- ObjectsProps(P).Centroid(2))^2);
            end
            %Calculate distance to closest neighbor
            Var.Measurements.(ObjName).(ImgName).ClosestNeighbor(O) = min(Distance([1:NumObjects]~=O));
            %Calculate number of neighbor within specified distance
            Var.Measurements.(ObjName).(ImgName).NumNeighbors(O) = length(find(Distance([1:NumObjects]~=O) <= Var.Analysis.NeighborhoodDistance));
        end
        
    end
    
    
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(floor(CallNum)) = toc;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = radtodeg(B)

A = B.*180./pi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
