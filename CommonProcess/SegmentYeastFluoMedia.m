function Var = SegmentYeastFluoMedia(Var, CallNum);
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
IncludeEdge = Var.SegmentationParameters(ParaNum).Touch;

SegmentMethod = Var.SegmentationParameters(ParaNum).SegmentMethod;
LocalMaximaType = Var.SegmentationParameters(ParaNum).LocalMaxType;
WatershedType = Var.SegmentationParameters(ParaNum).WatershedType;
ExcludeSize = Var.SegmentationParameters(ParaNum).ExcludeSize;

SegmentOptim = Var.SegmentationParameters(ParaNum).SegmentOptim;
HoleFillDisk = Var.SegmentationParameters(ParaNum).HoleFillDisk;
SegmentCleanDisk = Var.SegmentationParameters(ParaNum).SegmentCleanDisk;
Threshold = Var.SegmentationParameters(ParaNum).Threshold;



%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%


SegYFImg = Var.Analysis.SegYFImg{CallNum};

OrigImage = double(Var.Img.(SegYFImg));




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
    
    %Filter image
    sigma = 1;
    FiltLength = ceil(2*sigma);                                           % Determine filter size, min 3 pixels, max 61
    [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
    f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
    BlurredImage = conv2(OrigImage,f,'same');
    
    %Flatten the background of image based on expected cell size
    disk= round((MinDiam+MaxDiam)/2);
    mask        = strel('disk',disk);
    BlurredImage = imbothat(BlurredImage,mask);
    
    
    %Normalize image
    MinImg = min(BlurredImage(:));
    MaxImg = max(BlurredImage(:));
    NormImage = (BlurredImage-MinImg)/(MaxImg-MinImg);
%     NormImage = NormImage-mean(NormImage(:));
%     NormImage(NormImage<0) = 0;
    BlurredImage = NormImage;
    

    
    if strcmp(SegmentMethod, 'Edge')
        %Automatic edge detection
        [IntEdge, Threshold] = edge(BlurredImage,'sobel');
        if SegmentOptim ~= 1
            [IntEdge, Threshold] = edge(BlurredImage,'sobel', Threshold*SegmentOptim);
        end
    elseif strcmp(SegmentMethod, 'GreyThreshold')
        
        level = graythresh(BlurredImage);
        IntEdge = im2bw(BlurredImage,level*SegmentOptim);
    end
    
    
    %Clear object from borders
    IntEdge = imclearborder(IntEdge);
    
    % fill holes in image
    if HoleFillDisk ~= 0
        SE = strel('disk',HoleFillDisk);
        EdgeDilate = imdilate(IntEdge,SE);
        Edgefill = imfill(EdgeDilate,'holes');
        Objects = imerode(Edgefill,SE);
    end
    
    %remove small structures
    if SegmentCleanDisk ~= 0
        SE = strel('disk',SegmentCleanDisk);
        Objects = imopen(Objects,SE);
        Objects = imfill(Objects,'holes');
    end
    
    if strcmp(SegmentMethod, 'Edge')
        %%% Label the objects
        Objects = bwlabel(Objects);
        MeanBackground = mean(OrigImage(Objects == 0));
        NbObj = max(Objects(:));
        
        for Obj = 1:NbObj
            
            GrowObj = zeros(size(Objects));
            GrowObj(Objects == Obj) = 1;
            SE = strel('disk',15);
            GrowObj = imdilate(GrowObj,SE);
            MeanIntObj(Obj) = mean(OrigImage(Objects == Obj))-MeanBackground;
            MeanIntGrow(Obj) = mean(OrigImage(GrowObj == 1))-MeanBackground;
            if MeanIntObj(Obj) < 2*MeanIntGrow(Obj)
                Objects(Objects == Obj) = 0;
            end
        end
        if Debug == 1
            figure(103)
            plot(MeanIntObj, 'b')
            hold on
            plot(2.*MeanIntGrow, 'r')
            hold off
        end
    end
    
    
    if Debug == 1
        figure(100)
        subplot(2,2,1); imagesc(OrigImage), title('OrigImage') ; colormap(gray)
        subplot(2,2,2); imagesc(BlurredImage), title('BlurredImage') ; colormap(gray)
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
BlurredImage = conv2(padarray(NormImage, [FiltLength,FiltLength], 'replicate'),f,'same');
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
if ~isempty(LocalMaximaType)
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
    
    if Debug == 1
        figure(101)
        subplot(2,2,1); imagesc(NormImage), title('NormImage') ; colormap(gray)
        subplot(2,2,2); imagesc(Objects), title('Objects') ; colormap(gray)
        %         subplot(2,2,3); imagesc(Edgefill), title('Edgefill') ; colormap(gray)
        %         subplot(2,2,4); imagesc(Objects), title('Objects') ; colormap(gray)
        %         pause(1)
        
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

ObjProps = regionprops(Objects,'EquivDiameter', 'Area');
if strcmp(ExcludeSize,'Yes')
    DiamList =   [0;cat(1,ObjProps.EquivDiameter)];
    DiameterMap = DiamList(Objects+1);
    Objects(DiameterMap < MinDiam) = 0;
    Objects(DiameterMap > MaxDiam) = 0;
    
    %        AreaList =   [0;cat(1,ObjProps.Area)];
    %    AreaMap = AreaList(Objects+1);
    %    Objects(AreaMap < MinDiam) = 0;
    %     Objects(AreaMap > MaxDiam) = 0;
end


% %%% Get diameters of objects and calculate the interval
% %%% that contains 90% of the objects
% tmp = regionprops(Objects,'EquivDiameter');
% Diameters = [0;cat(1,tmp.EquivDiameter)];
% SortedDiameters = sort(Diameters);
% NbrInTails = max(round(0.05*length(Diameters)),1);
% Lower90Limit = SortedDiameters(NbrInTails);
% Upper90Limit = SortedDiameters(end-NbrInTails+1);
%
% %%% Locate objects with diameter outside the specified range
% tmp = Objects;
% if strcmp(ExcludeSize,'Yes')
%     %%% Create image with object intensity equal to the diameter
%     DiameterMap = Diameters(Objects+1);
%     %%% Remove objects that are too small
%     Objects(DiameterMap < MinDiam) = 0;
%     %%% Will be stored to the handles structure
%     SmallRemovedLabelMatrixImage = Objects;
%     %%% Remove objects that are too big
%     Objects(DiameterMap > MaxDiam) = 0;
% else
%     %%% Will be stored to the handles structure even if it's unedited.
%     SmallRemovedLabelMatrixImage = Objects;
% end
% %%% Store objects that fall outside diameter range for display
% DiameterExcludedObjects = tmp - Objects;

%%% Remove objects along the border of the image (depends on user input)
tmp = Objects;
if strcmp(IncludeEdge,'No')
    Objects = imclearborder(Objects);
end


%%% Relabel the objects
[Objects,NumOfObjects] = bwlabel(Objects > 0);

%% Display objects


if strcmp(Var.Figure.Display, 'on')

%ResultImg(:,:,1) = Objects > 0;
%ResultImg(:,:,2) = 0;
ResultImg(:,:,3) = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
ResultImg(:,:,1) = Objects > 0;


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
ObjProps = regionprops(Objects,'Centroid');

Var.Measurements.(SegYFImgOut).SegLabel = [1:length(ObjProps)];
if length(ObjProps) > 0
    for i = 1:length(ObjProps)
        Var.Measurements.(SegYFImgOut).CenterX(i,1) = ObjProps(i).Centroid(1);
        Var.Measurements.(SegYFImgOut).CenterY(i,1) = ObjProps(i).Centroid(2);
        Var.Measurements.(SegYFImgOut).SegLabel(i,1) = i;
    end
else
    error( 'No object found in in Image by Segmentation')
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;
