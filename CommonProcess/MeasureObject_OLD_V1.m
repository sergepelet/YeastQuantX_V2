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
       % O = O
        
        %Object Present in Frame
        %Size Features
        Var.Measurements.(ObjName).PixelList{1,Label(O)} = ObjectsProps(O).PixelIdxList;
        Var.Measurements.(ObjName).(ImgName).Area(Label(O)) = ObjectsProps(O).Area;
        Var.Measurements.(ObjName).(ImgName).Diameter(Label(O)) = ObjectsProps(O).EquivDiameter;
        Var.Measurements.(ObjName).(ImgName).CenterX(Label(O)) = ObjectsProps(O).Centroid(1);
        Var.Measurements.(ObjName).(ImgName).CenterY(Label(O)) = ObjectsProps(O).Centroid(2);
        
        Var.Measurements.(ObjName).(ImgName).CheckLabel(Label(O)) = Label(O);
        %Get label of parent object if it exists
        if isfield(Var.Measurements.(ObjName), 'ParentLabel')
            Var.Measurements.(ObjName).(ImgName).ParentLabel(Label(O)) = Var.Measurements.(ObjName).ParentLabel(O);
        end
        
        %Intensity features
        Var.Measurements.(ObjName).(ImgName).TotalIntensity(Label(O)) = sum(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).AverageIntensity(Label(O)) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MedianIntensity(Label(O)) = median(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).STDIntensity(Label(O)) = std(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MaxIntensity(Label(O)) = max(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).MinIntensity(Label(O)) = min(IntensityImage(ObjectsProps(O).PixelIdxList));
        Var.Measurements.(ObjName).(ImgName).STDIntensityNorm(Label(O)) = std(IntensityImage(ObjectsProps(O).PixelIdxList))./mean(IntensityImage(ObjectsProps(O).PixelIdxList));
        %HiPix calculation
        if ObjectsProps(O).Area <= NumHiPix
            Var.Measurements.(ObjName).(ImgName).LoPix(Label(O)) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            Var.Measurements.(ObjName).(ImgName).HiPix(Label(O)) = mean(IntensityImage(ObjectsProps(O).PixelIdxList));
            
            Var.Measurements.(ObjName).(ImgName).HiPixArea(Label(O)) = ObjectsProps(O).Area;
            
             Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(Label(O)) = NaN;
             Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(Label(O)) = NaN;
        else
            [SortedInt, SortedIndex] = sort(IntensityImage(ObjectsProps(O).PixelIdxList));
         
            Var.Measurements.(ObjName).(ImgName).LoPix(Label(O)) = mean(SortedInt(1:NumHiPix));
            Var.Measurements.(ObjName).(ImgName).HiPix(Label(O)) = mean(SortedInt(end-NumHiPix+1:end));
            
            %HiPix Area
            HiPixImg = zeros(size(ObjectImage));
            HiPixImg(ObjectsProps(O).PixelIdxList(SortedIndex(end-NumHiPix+1:end))) = 1;
            HiPixStat = regionprops(HiPixImg, 'ConvexArea');
            Var.Measurements.(ObjName).(ImgName).HiPixArea(Label(O)) = HiPixStat.ConvexArea;
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
                ConnectedHiPixStat = regionprops(LabelHP, IntensityImage, 'Area', 'MeanIntensity');  %'PixelIdxList',
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(Label(O)) = ConnectedHiPixStat.Area;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(Label(O)) = ConnectedHiPixStat.MeanIntensity;
                
            elseif NumArea == 0  %If no pixels are left after image opening set meas to NaN
                 Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(Label(O)) = NaN;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(Label(O)) = NaN;
            else   %If multiple objects are generated find the biggest one and measure its intensity
                CHPArea = 0;
                CHPInt = 0;
                for HPO = 1:NumArea
                    IndHPO = find(LabelHP == HPO);
                    if length(IndHPO) > CHPArea
                        CHPArea = length(IndHPO);
                        CHPInt = mean(IntensityImage(IndHPO));
                    end
                end
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPixArea(Label(O)) = CHPArea;
                Var.Measurements.(ObjName).(ImgName).ConnectedHiPix(Label(O)) = CHPInt;
            end
            
            
        end
        %Polarity Measurements using object shape
        Var.Measurements.(ObjName).(ImgName).Orientation(Label(O)) = 90-ObjectsProps(O).Orientation;
        Var.Measurements.(ObjName).(ImgName).MajorAxisLength(Label(O)) = ObjectsProps(O).MajorAxisLength;
        Var.Measurements.(ObjName).(ImgName).MinorAxisLength(Label(O)) = ObjectsProps(O).MinorAxisLength;
        Var.Measurements.(ObjName).(ImgName).Eccentricity(Label(O)) = ObjectsProps(O).Eccentricity;
        Var.Measurements.(ObjName).(ImgName).ConvexRatio(Label(O)) = ObjectsProps(O).Area/ObjectsProps(O).ConvexArea;
        
        %Polarity Measurements using weighted centroid
        Var.Measurements.(ObjName).(ImgName).WeightCenterX(Label(O)) = ObjectsProps(O).WeightedCentroid(1);
        Var.Measurements.(ObjName).(ImgName).WeightCenterY(Label(O)) = ObjectsProps(O).WeightedCentroid(2);
        WeightAxis = (ObjectsProps(O).WeightedCentroid(1)-ObjectsProps(O).Centroid(1))^2 + (ObjectsProps(O).WeightedCentroid(2)-ObjectsProps(O).Centroid(2))^2;
        Var.Measurements.(ObjName).(ImgName).WeightAxis(Label(O)) = WeightAxis;
        WeightAngle = atand(ObjectsProps(O).WeightedCentroid(1)-ObjectsProps(O).Centroid(1))/(ObjectsProps(O).WeightedCentroid(2)-ObjectsProps(O).Centroid(2));
        Var.Measurements.(ObjName).(ImgName).WeightAngle(Label(O)) = WeightAngle;
        
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
            Var.Measurements.(ObjName).(ImgName).Polarity(Label(O), 1:length(AngleList)-1 ) = PerimInt;

            [Var.Measurements.(ObjName).(ImgName).PolarityMaxInt(Label(O)), MaxInd] = max(PerimInt);
            Var.Measurements.(ObjName).(ImgName).PolarityMaxAngle(Label(O)) = Var.Analysis.AngleAxis(MaxInd);
            Var.Measurements.(ObjName).(ImgName).PolarityStdInt(Label(O)) = std(PerimInt);
            
            Var.Measurements.(ObjName).(ImgName).PerimSize(Label(O), 1:length(AngleList)-1 ) = PerimSize;
            
            Var.Measurements.(ObjName).(ImgName).Distance2Center(Label(O), 1:length(AngleList)-1 ) = Distance2Center;
            [Var.Measurements.(ObjName).(ImgName).D2C_Max(Label(O)), MaxInd] = max(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_MaxAngle(Label(O)) = Var.Analysis.AngleAxis(MaxInd);
            Var.Measurements.(ObjName).(ImgName).D2C_Mean(Label(O)) = mean(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_Std(Label(O)) = std(Distance2Center);
            Var.Measurements.(ObjName).(ImgName).D2C_StdNorm(Label(O)) = std(Distance2Center)/mean(Distance2Center);
            
            
            %% Use Circle object as center for polarity measurement
            CircleMethod = find(strcmp(Var.Analysis.SecMethod, 'Circle'));
            
            if ~isempty(CircleMethod)
            CircleObjName = Var.Analysis.SecObjOUT{CircleMethod};
                      
            %Get pixels relative to center of mass
            Xpix = ObjectsProps(O).PixelList(:,1) - Var.Measurements.(CircleObjName).CenterX(Label(O));
            Ypix = ObjectsProps(O).PixelList(:,2) - Var.Measurements.(CircleObjName).CenterY(Label(O));
            
            
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
            Var.Measurements.(ObjName).(ImgName).Polarity2CC(Label(O), 1:length(AngleList)-1 ) = CC_PerimInt;
            
            Var.Measurements.(ObjName).(ImgName).Distance2CC(Label(O), 1:length(AngleList)-1 ) = CC_Distance2Center;
             
            end
            
            
        end
    end
    
    
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(floor(CallNum)) = toc;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = radtodeg(B)

A = B.*180./pi;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
