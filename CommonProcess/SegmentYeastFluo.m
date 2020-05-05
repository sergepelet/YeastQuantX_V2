function Var = SegmentYeastFluo(Var, CallNum)

tic
Debug = 0;      %Set to 1 to display segmentation images and 0 not to

if nargin == 1
    CallNum = 1;
end


%Get Segmentation Parameters full name
SegParaFullName = [Var.Analysis.SegYFPara{CallNum}, '_', num2str(Var.Experiment.Objective),'x_bin', num2str(Var.Experiment.Bin)];

ParaNum = [];
%Identify correct segmentation paramters
for i = 1:length(Var.SegmentationParameters)
    if strcmp(SegParaFullName, Var.SegmentationParameters(i).FullName)
        ParaNum = i;
    end
end

if isempty(ParaNum)
    error('No proper  segmentation parameter found')
end

%load segmentation parameters
MinDiam = Var.SegmentationParameters(ParaNum).MinDiameter;
MaxDiam = Var.SegmentationParameters(ParaNum).MaxDiameter;
if isfield(Var.SegmentationParameters(ParaNum), 'DiameterType') && ~isempty(Var.SegmentationParameters(ParaNum).DiameterType)
    DiameterType = Var.SegmentationParameters(ParaNum).DiameterType;
else
    DiameterType = 'EquivDiameter';
end

MaxEccentricity = 0.9;
MinSolidity = 0.9;
IncludeEdge = Var.SegmentationParameters(ParaNum).Touch;
if isfield(Var.SegmentationParameters(ParaNum), 'RemoveHiIntPixels')
    RemovePix = Var.SegmentationParameters(ParaNum).RemoveHiIntPixels;
else
    RemovePix = 0;
end


SegmentMethod = Var.SegmentationParameters(ParaNum).SegmentMethod; % 
LocalMaximaType = Var.SegmentationParameters(ParaNum).LocalMaxType;
WatershedType = Var.SegmentationParameters(ParaNum).WatershedType;
ExcludeSize = Var.SegmentationParameters(ParaNum).ExcludeSize;
MedianFilter = Var.SegmentationParameters(ParaNum).MedianFilter;


SegmentOptim = Var.SegmentationParameters(ParaNum).SegmentOptim;
HoleFillDisk = Var.SegmentationParameters(ParaNum).HoleFillDisk;
SegmentCleanDisk = Var.SegmentationParameters(ParaNum).SegmentCleanDisk;
Threshold = Var.SegmentationParameters(ParaNum).Threshold;



%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%


SegYFImg = Var.Analysis.SegYFImg{CallNum};

OrigImage = double(Var.Img.(SegYFImg));
StartImg = OrigImage;
%if strcmp(MedianFilter, 'Yes')
OrigImage = medfilt2(OrigImage);
%end

HiThresh = prctile(OrigImage(:), 100-RemovePix);
OrigImage(OrigImage>HiThresh) = HiThresh;

figure(101)
imagesc(StartImg)
figure(102)
imagesc(OrigImage)

%%
% Define object in Image
if Threshold ~= 0
    %%% Apply a slight smoothing before thresholding to remove
    %%% 1-pixel objects and to smooth the edges of the objects.
    %%% Note that this smoothing is hard-coded, and not controlled
    %%% by the user.
    sigma = 1;
    FiltLength = ceil(2*sigma);                                           % Determine filter size, min 3 pixels, max 61
    [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
    f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
    BlurredImage = conv2(OrigImage,f,'same');                             % Blur original image
    Objects = BlurredImage > Threshold;                                   % Threshold image
    Threshold = mean(Threshold(:));                                       % Use average threshold downstreams
    Objects = imfill(double(Objects),'holes');                            % Fill holes
    
else
    %Normalize image
    MinImg = min(OrigImage(:));
    MaxImg = max(OrigImage(:));
    NormImage = (OrigImage-MinImg)/(MaxImg-MinImg);
    
    sigma = 1;
    FiltLength = ceil(2*sigma);                                           % Determine filter size, min 3 pixels, max 61
    [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
    f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
    BlurredImage = conv2(NormImage,f,'same');
    
    
    %Enhance image for objects of a given size range
    if strcmp(SegmentMethod, 'GreyThreshold_deltaBackground')
        %subtract background intensity from image
        %Open image with  disk ,uch larger than  objects
        BackgroundInt = imopen(BlurredImage,strel('disk',MaxDiam*10));
        %Gaussian Filter image
        BackgroundInt = imgaussfilt(BackgroundInt, MaxDiam*10/3);
        %Subtract Background intensity
        BlurredImage = BlurredImage-BackgroundInt;
        %Normalize image
        MinImg = min(BlurredImage(:));
        MaxImg = max(BlurredImage(:));
        BlurredImage = (BlurredImage-MinImg)/(MaxImg-MinImg);
        
    else
        disks=[MinDiam,MaxDiam];
        %disks=[50,100];
        for i=1:length(disks)
            mask        = strel('disk',disks(i));
            top         = imtophat(BlurredImage,mask);
            bot         = imbothat(BlurredImage,mask);
            BlurredImage    = imsubtract(imadd(BlurredImage,top), bot);
        end

    end
    
    if Debug
        figure(90)
        subplot(2,2,1); imagesc(OrigImage), title('OrigImage') ; colormap(gray)
        subplot(2,2,2); imagesc(NormImage), title('NormImage') ; colormap(gray)
        
        subplot(2,2,4); imagesc(BlurredImage), title('BlurredImage') ; colormap(gray)
        
        if strcmp(SegmentMethod, 'GreyThreshold_deltaBackground')
        subplot(2,2,3); imagesc(BackgroundInt), title('background') ; colormap(gray)
        end
        pause(1)
    end
    
    if strcmp(SegmentMethod, 'Edge')
        %Automatic edge detection
        [IntEdge, Threshold] = edge(BlurredImage,'sobel');
        if SegmentOptim ~= 1
            [IntEdge, Threshold] = edge(BlurredImage,'sobel', Threshold*SegmentOptim);
        end
        
        IntEdge(2,:) = 0;
        IntEdge(end-1,:) = 0;
        IntEdge(:,2) = 0;
        IntEdge(:,end-1) = 0;
        %            figure(99)
        %     imagesc(IntEdge)
        % IntEdge = imclearborder(IntEdge);
        
    elseif strcmp(SegmentMethod, 'GreyThreshold') || strcmp(SegmentMethod, 'GreyThreshold_deltaBackground')
        
%         %Segment image based on intensity
        level = graythresh(BlurredImage);
%         IntEdge = im2bw(BlurredImage,level*SegmentOptim);
        AllThresh = zeros(size(BlurredImage));
        %Test various thresholds
        for Fold = [1:0.1:SegmentOptim]
            if level*Fold<1
                %Calculate BW img
                testThresh = im2bw(BlurredImage,level*Fold);
                
                %Find object and their area
                testThresh = bwlabel(testThresh);
                ObjProps = regionprops(testThresh, DiameterType);
                
                %Filter objects based on their diameter size
                Diameters = [0;cat(1,ObjProps.(DiameterType))];
                DiameterMap = Diameters(testThresh+1);
                testThresh(DiameterMap < MinDiam) = 0;
                testThresh(DiameterMap > MaxDiam) = 0;
                %Add objects to final matrix
                AllThresh(testThresh>0) = AllThresh(testThresh>0) +1;
                
                if Debug >0
                    figure(round(Fold*1000))
                    imagesc(testThresh)
                    title(num2str(length(ObjProps)))
                end
                
            end
        end
        
        if Debug >0
            figure(120)
            imagesc(AllThresh)
        end
        %Transfer Objects to final segmentation image: Thersold set to 3...
        IntEdge = zeros(size(BlurredImage));
        IntEdge(AllThresh>3) = 1;
        
    end
    
    
    
    % fill holes in image
    if HoleFillDisk ~= 0
        SE = strel('disk',HoleFillDisk);
        EdgeDilate = imdilate(IntEdge,SE);
        %                figure(101)
        %     imagesc(EdgeDilate)
        Edgefill = imfill(EdgeDilate,'holes');
        Objects = imerode(Edgefill,SE);
    end
    %        figure(103)
    %     imagesc(Objects)
    %  drawnow
    %remove small structures
    if SegmentCleanDisk ~= 0
        SE = strel('disk',SegmentCleanDisk);
        Objects = imopen(Objects,SE);
        Objects = imfill(Objects,'holes');
    end
    
    %      figure(104)
    %     imagesc(Objects)
    
    if strcmp(SegmentMethod, 'Edge')
        %Remove objects whose intensity are not above the difference
        %between Backgroung and Objects
        MeanBackground = mean(OrigImage(EdgeDilate == 0));
        MeanObjects = mean(OrigImage(Objects > 0));
        Delta = MeanObjects-MeanBackground;
        Objects(OrigImage<MeanBackground+0.8*Delta) = 0;
        
        %       figure(105)
        %     imagesc(Objects); title('remove Delta')
        
        
        %Get Eccentricity and Solidity from objects
        LabelObj = bwlabel(Objects);
        Props = regionprops(LabelObj, 'Eccentricity', 'Solidity'); %, 'MajorAxisLength', 'MeanIntensity');
        %Create Soliditiy Map
        Solidity = [0;cat(1,Props.Solidity)];
        SolidityyMap = Solidity(LabelObj+1);
        %remove objects below Minimal Solidity
        Objects(SolidityyMap < MinSolidity) = 0;
        
        %Create Eccentricity Map
        Eccentricity = [0;cat(1,Props.Eccentricity)];
        EccentricityMap = Eccentricity(LabelObj+1);
        %remove objects with large eccentricitiy
        Objects(EccentricityMap >MaxEccentricity) = 0;
        figure(107)
        imagesc(SolidityyMap); title('SolidityyMap')
        %
        %             figure(106)
        %     imagesc(Objects); title('remove Eccentricity')
        % MeanInt = [0;cat(1,Props.MeanIntensity)];
        % MeanIntMap = MeanInt(LabelObj+1);
        %MajorAxis = [0;cat(1,Props.MajorAxisLength)];
        % MajorAxisMap = MajorAxis(LabelObj+1);
        
        %           figure(107)
        %         imagesc(EccentricityMap); title('EccentricityMap')
        %
        %         figure(107)
        %         imagesc(MajorAxisMap); title('MajorAxisMap')
        %         figure(108)
        %         imagesc(MeanIntMap); title('MeanIntMap')
        %
    end
    
    %       figure(105)
    %     imagesc(Objects)
    
    %     if strcmp(SegmentMethod, 'Edge')
    %         %%% Label the objects
    %         Objects = bwlabel(Objects);
    %         MeanBackground = mean(OrigImage(Objects == 0));
    %         NbObj = max(Objects(:));
    %
    %             for Obj = 1:NbObj
    %
    %                 GrowObj = zeros(size(Objects));
    %                 GrowObj(Objects == Obj) = 1;
    %                 SE = strel('disk',15);
    %                 GrowObj = imdilate(GrowObj,SE);
    %                 MeanIntObj(Obj) = mean(OrigImage(Objects == Obj))-MeanBackground;
    %                 MeanIntGrow(Obj) = mean(OrigImage(GrowObj == 1))-MeanBackground;
    %                 if MeanIntObj(Obj) < 2*MeanIntGrow(Obj)
    %                     Objects(Objects == Obj) = 0;
    %                 end
    %             end
    %             if Debug == 1
    %             figure(103)
    %             plot(MeanIntObj, 'b')
    %             hold on
    %             plot(2.*MeanIntGrow, 'r')
    %             hold off
    %             end
    %     end
    %
    %
    %     %Make sure that edge have not been split in many objects if it is the
    %     %case remove the low intensity ones
    %     if strcmp(SegmentMethod, 'Edge')
    %         %%% Label the objects
    %         Objects = bwlabel(Objects);
    %         Edgefill = bwlabel(Edgefill);
    %         FillProps = regionprops(Edgefill, 'PixelIdxList');
    %         NbObj = length(FillProps);
    %         if NbObj>0
    %             for i = 1:NbObj
    %                 ObjLabel = unique(Objects(FillProps(i).PixelIdxList));
    %                 NumObj = find(ObjLabel>0);
    %                 if NumObj >1
    %                     MeanInt = zeros(length(ObjLabel),1);
    %                     for j = 1:length(ObjLabel)
    %                         if ObjLabel > 0
    %
    %                             MeanInt(j) = mean(OrigImage(Objects == ObjLabel(j)));
    %                         end
    %
    %                     end
    %                     [SortInt, SortIndex] = sort(MeanInt);
    %                     for j = 1:length(ObjLabel)-1
    %                     	Objects(Objects == ObjLabel(SortIndex(j))) = 0;
    %                     end
    %                 end
    %             end
    % %             MeanAll = mean(MeanObjInt);
    % %             MeanBackground = mean(OrigImage(Objects == 0));
    % %
    % %             LowIntObj = find((MeanObjInt-MeanBackground) < (MeanAll-MeanBackground)*0.15)
    % %             if length(LowIntObj)>0
    % %                 for j = 1:length(LowIntObj)
    % %                     Objects(Objects == LowIntObj(j)) = 0;
    % %                 end
    % %                 MeanAll = MeanAll
    % %                 MeanBackground = MeanBackground
    % %                 MeanObjInt = MeanObjInt
    % %             end
    %         end
    %     end
    
    
    if Debug
        figure(100)
        subplot(2,2,1); imagesc(NormImage), title('NormImage') ; colormap(gray)
        subplot(2,2,2); imagesc(IntEdge), title('Edges') ; colormap(gray)
        subplot(2,2,3); imagesc(Edgefill), title('Edgefill') ; colormap(gray)
        subplot(2,2,4); imagesc(Objects), title('Objects') ; colormap(gray)
        pause(1)
        
    end
end



%Segment touching objects

sigma = MinDiam/3.5;
FiltLength = min(30,max(1,ceil(2*sigma)));                            % Determine filter size, min 3 pixels, max 61
[x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
%%% The original image is blurred. Prior to this blurring, the
%%% image is padded with values at the edges so that the values
%%% around the edge of the image are not artificially low.  After
%%% blurring, these extra padded rows and columns are removed.
%BlurredImage = conv2(padarray(GFPImage, [FiltLength,FiltLength], 'replicate'),f,'same');
BlurredImage = conv2(padarray(OrigImage, [FiltLength,FiltLength], 'replicate'),f,'same');
BlurredImage = BlurredImage(FiltLength+1:end-FiltLength,FiltLength+1:end-FiltLength);
%%% Get local maxima, where the definition of local depends on the
%%% user-provided object size. This will (usually) be done in a
%%% lower-resolution image for speed. The ordfilt2() function is
%%% very slow for large images containing large objects.
%%% Therefore, image is resized to a size where the smallest
%%% objects are about 10 pixels wide. Local maxima within a radius
%%% of 5-6 pixels are then extracted. It might be necessary to
%%% tune this parameter. The MaximaSuppressionSize must be an
%%% integer.  The MaximaSuppressionSize should be equal to the
%%% minimum acceptable radius if the objects are perfectly
%%% circular with local maxima in the center. In practice, the
%%% MinDiameter is divided by 1.5 because this allows the local
%%% maxima to be shifted somewhat from the center of the object.

MaximaSuppressionSize = round(MinDiam/1.5);
MaximaMask = getnhood(strel('disk', MaximaSuppressionSize));
if ~strcmp(LocalMaximaType,'None')
    if strcmp(LocalMaximaType,'Intensity')
        
        MaximaMask = getnhood(strel('disk', min(50,max(1,floor(MinDiam/1.5)))));
        %Initialize MaximaImage
        MaximaImage = BlurredImage;
        %Save only local maxima
        MaximaImage(BlurredImage < ...
            ordfilt2(BlurredImage,sum(MaximaMask(:)),MaximaMask)) = 0;
        %Remove dim maxima
        MaximaImage = MaximaImage > Threshold;
        
    elseif strcmp(LocalMaximaType,'Shape')
        %%% Calculate distance transform
        DistanceTransformedImage = bwdist(~Objects);
        %%% Add some noise to get distinct maxima
        DistanceTransformedImage = DistanceTransformedImage + ...
            0.001*rand(size(DistanceTransformedImage));
        ImageResizeFactor = 1;
        ResizedDistanceTransformedImage = imresize(DistanceTransformedImage,ImageResizeFactor,'bilinear');
        %%% Initialize MaximaImage
        MaximaImage = ones(size(ResizedDistanceTransformedImage));
        %%% Set all pixels that are not local maxima to zero
        MaximaImage(ResizedDistanceTransformedImage < ...
            ordfilt2(ResizedDistanceTransformedImage,sum(MaximaMask(:)),MaximaMask)) = 0;
        %%% Restore image size
        MaximaImage = imresize(MaximaImage,size(Objects),'bilinear');
        %%% We are only interested in maxima within thresholded objects
        MaximaImage(~Objects) = 0;
        %%% Shrink to points (needed because of the resizing)
        MaximaImage = bwmorph(MaximaImage,'shrink',inf);
    end
    
    %%% Overlay the maxima on either the original image or a distance
    %%% transformed image. The watershed is currently done on
    %%% non-smoothed versions of these image. We may want to try to do
    %%% the watershed in the slightly smoothed image.
    if strcmp(WatershedType,'Intensity')
        %%% Overlays the objects markers (maxima) on the inverted original image so
        %%% there are black dots on top of each dark object on a white background.
        Overlaid = imimposemin(1 - OrigImage,MaximaImage);
    elseif strcmp(WatershedType,'Membrane')
        
        MedSize= 5;
        FiltImg = medfilt2(OrigImage, [MedSize MedSize]);
        Overlaid = imimposemin(FiltImg,MaximaImage);
        %Overlaid = FiltImg;
        Overlaid(Objects==0)=Inf;
        
    elseif strcmp(WatershedType,'Distance')
        %%% Overlays the object markers (maxima) on the inverted DistanceTransformedImage so
        %%% there are black dots on top of each dark object on a white background.
        %%% We may have to calculate the distance transform:
        if ~exist('DistanceTransformedImage','var')
            DistanceTransformedImage = bwdist(~Objects);
        end
        Overlaid = imimposemin(-DistanceTransformedImage,MaximaImage);
    end
    
    %%% Calculate the watershed transform and cut objects along the boundaries
    WatershedBoundaries = watershed(Overlaid) > 0;
    Objects = Objects.*WatershedBoundaries;
    
    
    if strcmp(WatershedType,'Membrane')
        
        SE = strel('disk',10);
        MaximaImage = imdilate(MaximaImage,SE);
        
        CombineImg = Objects;
        CombineImg(MaximaImage) = 2;
        
        MainObj = imextendedmax(CombineImg,1);
        
        
        MainDist = bwdist(MainObj);
        WaterMain = watershed(MainDist);
        
        
        [LableObj, NumObj] = bwlabel(MainObj);
        SE = strel('disk',SegmentCleanDisk);
        NewMainObj = zeros(size(MainObj));
        for O = 1:NumObj
            SingleObj = zeros(size(MainObj));
            SingleObj(LableObj== O) = 1;
            
            IndSingle = find(SingleObj == 1);
            WaterVal = unique(WaterMain(IndSingle));
            if length(WaterVal) > 1
                AllVal = WaterMain(IndSingle);
                NumInd = [];
                for i = 1:length(WaterVal)
                    NumInd(i) = length(find(AllVal == WaterVal(i)));
                end
                [SortNumInd, IndSort] = sort(NumInd, 'descend');
                WaterVal = WaterVal(IndSort(1));
            end
            
            
            SingleObj = imdilate(SingleObj,SE);
            SingleObj = imfill(SingleObj, 'holes');
            SingleObj = imerode(SingleObj,SE);
            
            
            SE = strel('disk',4);
            SingleObj = imdilate(SingleObj,SE);
            SE = strel('disk',SegmentCleanDisk);
            SingleObj = imopen(SingleObj,SE);
            
            SingleObj(WaterMain ~= WaterVal) = 0;
            
            [LableSingleObj, NumSingleObj] =bwlabel(SingleObj);
            if NumSingleObj ~= 1
                NumSingleObj = NumSingleObj
            end
            NewMainObj(SingleObj == 1) = O;
        end
        
        
        Objects = NewMainObj;
         if Debug
        figure(101)
        Perim = bwperim(NewMainObj);
        FiltImg(Perim == 1) = 500;
        subplot(2,2,1); imagesc(CombineImg), title('CombineImg')
        subplot(2,2,4); imagesc(NewMainObj), title('NewMainObj')
        subplot(2,2,3); imagesc(FiltImg), title('FiltImg')
        subplot(2,2,2); imagesc(MainObj), title('MainObj')
        pause(1)
        
    end
    end
    
   
    
    
    %%% Label the objects
    Objects = bwlabel(Objects);
    
    %%% Remove objects with no marker in them (this happens occasionally)
    %%% This is a very fast way to get pixel indexes for the objects
    tmp = regionprops(Objects,'PixelIdxList');
    for k = 1:length(tmp)
        %%% If there is no maxima in these pixels, exclude object
        if sum(MaximaImage(tmp(k).PixelIdxList)) == 0
            Objects(tmp(k).PixelIdxList) = 0;
        end
    end
    
end
drawnow

%%% Label the objects
Objects = bwlabel(Objects);

AllObjects = Objects;

% %%% Merge small objects
% if strcmp(MergeChoice,'Yes')
%     Objects = MergeObjects(Objects,OrigImage,[MinDiameter MaxDiameter]);
% end


%%% Get diameters of objects and calculate the interval
%%% that contains 90% of the objects
tmp = regionprops(Objects,DiameterType);
Diameters = [0;cat(1,tmp.(DiameterType))];
SortedDiameters = sort(Diameters);
NbrInTails = max(round(0.05*length(Diameters)),1);
Lower90Limit = SortedDiameters(NbrInTails);
Upper90Limit = SortedDiameters(end-NbrInTails+1);

%%% Locate objects with diameter outside the specified range
tmp = Objects;
if Debug
    figure(110)
    imagesc(Objects); title('All segmented Objects')
end
    
if strcmp(ExcludeSize,'Yes')
    %%% Create image with object intensity equal to the diameter
    DiameterMap = Diameters(Objects+1);
%     figure(100)
%     imagesc(DiameterMap)
    %%% Remove objects that are too small
    Objects(DiameterMap < MinDiam) = 0;
    %%% Will be stored to the handles structure
    SmallRemovedLabelMatrixImage = Objects;
    %%% Remove objects that are too big
    Objects(DiameterMap > MaxDiam) = 0;
else
    %%% Will be stored to the handles structure even if it's unedited.
    SmallRemovedLabelMatrixImage = Objects;
end
%%% Store objects that fall outside diameter range for display
DiameterExcludedObjects = tmp - Objects;

%%% Remove objects along the border of the image (depends on user input)
tmp = Objects;
if strcmp(IncludeEdge,'No')
    Objects = imclearborder(Objects);
end


%%% Relabel the objects
[Objects,NumOfObjects] = bwlabel(Objects > 0);




%ResultImg(:,:,1) = Objects > 0;

ResultImg(:,:,2) = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
ResultImg(:,:,1) = Objects > 0;
ResultImg(:,:,3) = 0;


% % %%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Bypass everything:
% % %%% Relabel the objects
% % [Objects,NumOfObjects] = bwlabel(Var.Img.GroupSegYeast > 0);
% %
%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow



if strcmp(Var.Figure.Display, 'on')
    FigNum = find(strcmp(Var.Figure.List, 'SegYF'));
    figure(FigNum(CallNum))
    
    %%% A subplot of the figure window is set to display the original image.
    subplot(2,2,1); imagesc(OrigImage); title('Input Image');
    %%% A subplot of the figure window is set to display the colored label
    %%% matrix image.
    subplot(2,2,2); imagesc(AllObjects); title('All Objects found');
    
    %%% A subplot of the figure window is set to display the inverted original
    %%% image with watershed lines drawn to divide up clusters of objects.
    subplot(2,2,4); imagesc(Objects); title([num2str(NumOfObjects), ' Retained Objetcs']);
    
    %%% A subplot of the figure window is set to display the inverted original
    %%% image with watershed lines drawn to divide up clusters of objects.
    subplot(2,2,3); image(ResultImg); title([num2str(NumOfObjects), 'Objetcs on Image']);
    
    %pause(5)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow



%%% Saves the segmented image

SegYFImgOut = Var.Analysis.SegYFImgOut{CallNum};

Var.Img.(SegYFImgOut) = Objects;
%Save X-Y center position of objects
ObjProps = regionprops(Objects,'Centroid','EquivDiameter');
if length(ObjProps) > 0
    for i = 1:length(ObjProps)
        Var.Measurements.(SegYFImgOut).CenterX(i,1) = ObjProps(i).Centroid(1);
        Var.Measurements.(SegYFImgOut).CenterY(i,1) = ObjProps(i).Centroid(2);
        Var.Measurements.(SegYFImgOut).EquivDiameter(i,1) = ObjProps(i).EquivDiameter;
        Var.Measurements.(SegYFImgOut).SegLabel(i,1) = i;
    end
else
    error( 'No object found in in Image by Segmentation')
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;

