function Var = SegmentYeastPhase(Var, CallNum)
tic
Debug = 0;     %Set to 1 to display segmentation images and 0 not to
DebugIter = 1;


%Get Segmentation Parameters full name
SegParaFullName = [Var.Analysis.SegPhasePara{CallNum}, '_', num2str(Var.Experiment.Objective),'x_bin', num2str(Var.Experiment.Bin)];

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


%min and Max diameter of objects
MinDiam = Var.SegmentationParameters(ParaNum).MinDiameter;



%Smallest area for objects
ThreshArea = Var.SegmentationParameters(ParaNum).MinArea;
%Max area for objects
MaxArea = Var.SegmentationParameters(ParaNum).MaxArea;


%Values used for DIC image optimization
SmallValue1 = Var.SegmentationParameters(ParaNum).SmallBlur;
LargeValue1 = Var.SegmentationParameters(ParaNum).LargeBlur;

%Values for dilation of borders
SmallDisk =  Var.SegmentationParameters(ParaNum).SmallDilation;
MedDisk = Var.SegmentationParameters(ParaNum).LargeDilation;

% Threshold factor
MultThresh =  1.0; %2.0; %
%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%

%Transfer DIC image to OrigImage variable
OrigImage = double(Var.Img.(Var.Analysis.SegPhiImg{CallNum}));

%Check if second DIC image at a different focus
if isfield(Var.Analysis, 'SegPhiSecImg') && ~isempty(Var.Analysis.SegPhiSecImg);
    SecImage = double(Var.Img.(Var.Analysis.SegPhiSecImg{CallNum}));
    SecondZ = 1;
else
    SecondZ = 0;
end

% %Read previous segmentation from Var if not the first frame
% if ~(strcmp(Var.Analysis.FrameAnalyzed, 'First') || strcmp(Var.Analysis.FrameAnalyzed, 'Single'))
%     PrevSeg = Var.Img.(Var.Analysis.SegPhiImgOut{CallNum});
% end
% 


%% Enhance image  for border detection
%Remove Hi or Low Intensity Pix
DeviationFold = 2;
MeanInt = mean(OrigImage(:));
StdInt = std(OrigImage(:));
LowPix = find(OrigImage < MeanInt - DeviationFold*StdInt);
OrigImage(LowPix) = MeanInt - DeviationFold*StdInt;
HiPix = find(OrigImage > MeanInt + DeviationFold*StdInt);
OrigImage(HiPix) = MeanInt + DeviationFold*StdInt;


%Blurr image
sigma = 1;
FiltLength = ceil(2*sigma);                                           % Determine filter size, min 3 pixels, max 61
[x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
BlurredImage = conv2(OrigImage -mean(OrigImage(:)),f,'same');                             % Blur original image

%Enhance image for objects of a given size range

disks=[SmallValue1 LargeValue1];
for i=1:length(disks)
    mask        = strel('disk',disks(i));
    top         = imtophat(BlurredImage,mask);
    bot         = imbothat(BlurredImage,mask);
    EnhancedInvertedImage    = imsubtract(imadd(BlurredImage,top), bot);
    drawnow
end

%Fill Image
EnhancedZ1 = imfill(abs(EnhancedInvertedImage));

%Automatic edge detection
[BWLowEdgeZ1, ThresholdZ1] = edge(EnhancedZ1,'sobel');
BWLowEdge = BWLowEdgeZ1;
%Precess the second image in the same fashion
if SecondZ
    %Threshold Hi or Low Intensity Pix
    MeanInt = mean(SecImage(:));
    StdInt = std(SecImage(:));
    LowPix = find(SecImage < MeanInt - DeviationFold*StdInt);
    SecImage(LowPix) = MeanInt - DeviationFold*StdInt;
    HiPix = find(SecImage > MeanInt + DeviationFold*StdInt);
    SecImage(HiPix) = MeanInt + DeviationFold*StdInt;
    
    %Blur
    sigma = 1;
    FiltLength = ceil(2*sigma);                                           % Determine filter size, min 3 pixels, max 61
    [x,y] = meshgrid(-FiltLength:FiltLength,-FiltLength:FiltLength);      % Filter kernel grid
    f = exp(-(x.^2+y.^2)/(2*sigma^2));f = f/sum(f(:));                    % Gaussian filter kernel
    BlurredImage = conv2(SecImage -mean(SecImage(:)),f,'same');
    disks=[SmallValue1 LargeValue1];
    for i=1:length(disks)
        mask        = strel('disk',disks(i));
        top         = imtophat(BlurredImage,mask);
        bot         = imbothat(BlurredImage,mask);
        EnhancedInvertedImage    = imsubtract(imadd(BlurredImage,top), bot);
        drawnow
    end
    %Fill Image
    EnhancedZ2 = imfill(abs(EnhancedInvertedImage));
    
    %Automatic edge detection
    [BWLowEdgeZ2, ThresholdZ2] = edge(EnhancedZ2,'sobel');
    
    %add eges from second Z image to first image
    BWLowEdge(BWLowEdgeZ2==1) = 1;
end
if Debug
    figure(200)
    imagesc(BWLowEdge)
end
%Dilation and erosion using a Small disk on low edge to maximize cell area in the
%segmentation
SE = strel('disk',SmallDisk);
BWDilate = imdilate(BWLowEdge,SE);
BWFill = imfill(BWDilate,'holes');
BWLowGroupCell = imerode(BWFill,SE);

%Clean  up cell contour using medium disk
SE = strel('disk',MedDisk);
BWLowGroupCell = imopen(BWLowGroupCell, SE);
if Debug
    figure(201)
    imagesc(BWLowGroupCell)
end
%% Define High edge image by adjusting the threshold for the segmentation

%Loop until the segmentation threshold allows a good match between present
%image and previous segmentation
MatchPrevSeg = 0;
while ~MatchPrevSeg
    %
    [BWEdgeZ1, HiThresholdZ1] = edge(EnhancedZ1,'sobel', ThresholdZ1*MultThresh);
    
    
    %Define Edge image used for filling
    %Use a combination of the two edge image found
    BWEdge = BWEdgeZ1;
    %Use a combination of the two edges image if 2 focal plane
    if SecondZ
        [BWEdgeZ2, HiThresholdZ2] = edge(EnhancedZ2,'sobel', ThresholdZ2*MultThresh);
        BWEdge(BWEdgeZ2==1) = 1;
    end
    
    %Remove edges toughing the border
    BWEdge = imclearborder(BWEdge);
    
    %Dilation and erosion using a small disk to minimize artifact in the
    %segmentation
    SE = strel('disk',SmallDisk);
    BWDilate = imdilate(BWEdge,SE);
    BWFill = imfill(BWDilate,'holes');
    BWGroupCell = imerode(BWFill,SE);
    BWGroupCell = imopen(BWGroupCell, SE);
    if Debug
        figure(202)
        imagesc(BWGroupCell)
    end
    
    %Label Area found with high edge and define properites
    [LabelGroupCell,NumGroupCell] = bwlabel(BWGroupCell);
    Groupprops = regionprops(LabelGroupCell,'Area','MajorAxisLength', 'MinorAxisLength', 'PixelIdxList');
    
    
    %check that aea and diameter is not too small to be a cell
    for m = 1:NumGroupCell
        if Groupprops(m).MajorAxisLength < MinDiam  || Groupprops(m).Area < ThreshArea
            %             MajorAxis = Groupprops(m).MajorAxisLength
            %             Area = Groupprops(m).Area
            BWGroupCell(Groupprops(m).PixelIdxList) = 0;
        end
    end
    RawCombinedGroups = BWGroupCell;
    
    %Combine two group images
    CombinedGroups = double(BWLowGroupCell);
    CombinedGroups(BWGroupCell>0) = 2;
    
    %Keep cells  shape from low edge, but only if they overlap with cells found
    %with the high threshold
    CombinedGroups = imextendedmax(CombinedGroups,1);
    
    [LabelGroupCell,NumGroupCell] = bwlabel(CombinedGroups);
    Groupprops = regionprops(LabelGroupCell,'Area','MajorAxisLength', 'MinorAxisLength', 'PixelIdxList');
    
    for m = 1:NumGroupCell
        if Groupprops(m).Area > MaxArea
            %  LargeArea = Groupprops(m).Area
            CombinedGroups(Groupprops(m).PixelIdxList) = 0;
        end
    end
    
    
    %If there is a previous segmentation image
    if exist('PrevSeg')
        %Combine two group images
        CombinedSegImg = double(CombinedGroups);
        CombinedSegImg(PrevSeg>0) = CombinedSegImg(PrevSeg>0)+1;
        
        %Get cell that overlap each other
        OverlapSegImg = imextendedmax(CombinedSegImg,1);
        
        %Count cells in images
        [LabelCombinedSegImg,NumOverlapCell] = bwlabel(OverlapSegImg);
        [LabelPrevSeg,NumPrevCell] = bwlabel(PrevSeg);
        
        if Debug
            figure(100+DebugIter)
            DebugIter = DebugIter+1;
            subplot(2,2,1); imagesc(CombinedGroups), title('CombinedGroups')
            subplot(2,2,2); imagesc(CombinedSegImg), title('CombinedSegImg')
            subplot(2,2,3); imagesc(PrevSeg), title('PrevSeg')
            subplot(2,2,4); imagesc(OverlapSegImg), title('OverlapSegImg')
        end
        
        
        %                     NumGroupCell = NumGroupCell
        %             NumPrevCell = NumPrevCell
        %             NumOverlapCell = NumOverlapCell
        %If number of cell in overlap image is equal to present and past cell
        %numbers segmentation is good
        if NumGroupCell == NumOverlapCell && NumPrevCell == NumOverlapCell
            MatchPrevSeg = 1;
            %If there are less cell in the overlap image than in the previous
            %segmentation image try to lower the threshold for segmentation
        elseif NumOverlapCell < NumPrevCell
            
            MultThresh = MultThresh-0.25;
            if MultThresh < 1.25
                %If a very low threshold cannot do the job then forget
                %about it
                MatchPrevSeg = 1;
            else
                MatchPrevSeg = 0;
            end
        else
            
            MatchPrevSeg = 1;
        end
        
        
    else
        MatchPrevSeg = 1;
    end
end




if Debug
    EnhancedAll = EnhancedZ1 + EnhancedZ2;
    EnhancedRGB(:,:,1) = EnhancedZ1./max(EnhancedZ1(:));
    EnhancedRGB(:,:,2) = EnhancedZ2./max(EnhancedZ2(:));
    EnhancedRGB(:,:,3) = 0;
    figure(100+DebugIter)
    DebugIter = DebugIter+1;
    subplot(2,2,1); imagesc(EnhancedZ1), title('EnhancedZ1')
    subplot(2,2,2); imagesc(EnhancedZ2), title('EnhancedZ2')
    subplot(2,2,3); imagesc(EnhancedRGB), title('EnhancedRGB')
    subplot(2,2,4); imagesc(EnhancedAll), title('EnhancedAll')
    
    BWEdgeRGB(:,:,1) = BWEdgeZ1;
    BWEdgeRGB(:,:,2) = BWEdgeZ2;
    BWEdgeRGB(:,:,3) = 0;
    figure(100+DebugIter)
    DebugIter = DebugIter+1;
    subplot(2,2,1); imagesc(BWEdgeZ1), title('BWEdgeZ1')
    subplot(2,2,2); imagesc(BWEdgeZ2), title('BWEdgeZ2')
    %subplot(2,2,3); imagesc(BWEdgeZA), title('BWEdgeZA')
    subplot(2,2,4); imagesc(BWEdgeRGB), title('BWEdgeRGB')
end

if Debug
    GroupR = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    GroupR(BWGroupCell) = 0.5;
    GroupG = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    GroupG(BWLowGroupCell) = 0.5;
    GroupRGB(:,:,1) = GroupR;
    GroupRGB(:,:,2) = GroupG;
    GroupRGB(:,:,3) = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    
    figure(100+DebugIter)
    DebugIter = DebugIter+1;
    subplot(2,2,1); imagesc(BWLowGroupCell), title('BWLowGroupCell')
    subplot(2,2,2); imagesc(RawCombinedGroups), title('RawCombinedGroups')
    subplot(2,2,3); image(GroupRGB), title('GroupRGB')
    subplot(2,2,4); imagesc(CombinedGroups), title('CombinedGroups')
    
end

%Clean  up cell contour using medium disk
SE = strel('disk',MedDisk);
CombinedGroups = imopen(CombinedGroups, SE);
if Debug
    figure(203)
    imagesc(CombinedGroups)
end
%Label the cell
[LabelGroupCell,NumGroupCell] = bwlabel(CombinedGroups);
Groupprops = regionprops(LabelGroupCell,'BoundingBox','ConvexArea','ConvexImage', 'Area');

%Remove from edge image the elements that don't belong to a group of cell
BWEdge(CombinedGroups==0) = 0;

%%
% %Border parameter sets the number of pixels added on each sides of the
% %bounding box defined around the large areas
% Border = 5;
%
% FinalGroupCell = zeros(size(BWGroupCell));
% if NumGroupCell > 0
%     for n = 1:NumGroupCell
%         Create image containing only one large area
%         use the bounding box to define the corners of the small image
%         ULCorner = round(Groupprops(n).BoundingBox(1:2));
%         Width = Groupprops(n).BoundingBox(3:4);
%         Corners = [ULCorner(2),ULCorner(2)+Width(2),ULCorner(1),ULCorner(1)+Width(1) ];
%         Check if corners are within the size of the full image
%         if ULCorner(1) + Width(1) > size(BWGroupCell,2)
%             Corners(4) = size(BWGroupCell,2);
%         end
%         if ULCorner(2) + Width(2) > size(BWGroupCell,1)
%             Corners(2) = size(BWGroupCell,1);
%         end
%
%         Cut full images to keep only the small region that contain the
%         area of interest
%         CutGroup = CutWithBorder(LabelGroupCell, Corners, Border);
%         CutGroup(CutGroup~= n) = 0;
%
%         CutHull = zeros(size(CutGroup));
%                 SC = size(CutHull(Border:end-Border,Border:end-Border))
%                 SH = size(Groupprops(n).ConvexImage)
%
%         CutHull(Border+1:end-Border-1,Border+1:end-Border-1) = Groupprops(n).ConvexImage;
%
%         Check if convex area is much larger than groupCell.
%
%         if Groupprops(n).Area < 0.5* Groupprops(n).ConvexArea
%
%             Do the fill process with the low edge image
%             CutEdge = CutWithBorder(BWLowEdge, Corners, Border);
%             SE = strel('disk',SmallDisk);
%             BWDilate = imdilate(CutEdge,SE);
%             BWCutCell = imfill(CutEdge,'holes');
%             BWCutCell = imerode(BWFill,SE);
%             if Debug
%                 SizeComp = BWCutCell;
%                 SizeComp(CutGroup>0) = SizeComp(CutGroup>0) + 1;
%                 SizeComp(CutHull>0) = SizeComp(CutHull>0) + 2;
%                 figure(200 +n); imagesc(SizeComp);title(['Area:', num2str(Groupprops(n).Area), 'Hull:', num2str(Groupprops(n).ConvexArea)]);
%                 pause(4)
%             end
%
%             CutGroup = BWCutCell.*n;
%
%         else
%
%             if Debug
%                 SizeComp = CutGroup;
%
%                 SizeComp(CutHull>0) = SizeComp(CutHull>0) + 2;
%                 figure(200 +n); imagesc(SizeComp); title(['Area:', num2str(Groupprops(n).Area), 'Hull:', num2str(Groupprops(n).ConvexArea)]);
%                 pause(4)
%             end
%
%
%         end
%
%         Perform image closeing on cut image to smooth it
%         SE = strel('disk',2);
%         CutGroup = imclose(CutGroup,SE);
%
%         Transfer group cell in Final image
%         FinalGroupCell(Corners(1):Corners(2),Corners(3):Corners(4)) = FinalGroupCell(Corners(1):Corners(2),Corners(3):Corners(4)) + CutGroup(Border+1:end-Border,Border+1:end-Border);
%
%     end
% end

% if Debug
%     figure(107)
%     subplot(2,2,1); imagesc(BWGroupCell), title('BWGroupCell')
%     subplot(2,2,2); imagesc(LabelGroupCell), title('LabelGroupCell')
%     subplot(2,2,3); imagesc(FinalGroupCell), title('FinalGroupCell')
%     subplot(2,2,4); imagesc(FinalGroupCell), title('FinalGroupCell')
% end


%%


%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%
drawnow



if strcmp(Var.Figure.Display, 'on')
    
    GroupR = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    GroupR(bwperim(LabelGroupCell)) = 1;
    GroupG = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    %GroupG(BWLowGroupCell) = 0.5;
    GroupRGB(:,:,1) = GroupR;
    GroupRGB(:,:,2) = GroupG;
    GroupRGB(:,:,3) = (OrigImage-min(OrigImage(:)))./(max(OrigImage(:))-min(OrigImage(:)));
    
    FigNum = find(strcmp(Var.Figure.List, 'SegPhi'));
    figure(FigNum(CallNum))
    subplot(2,2,1); imagesc(OrigImage); title('Input Image');
    subplot(2,2,2); imagesc(BWEdge); title('BWEdge');
    subplot(2,2,3); imagesc(LabelGroupCell); title('LabelGroupCell');
    subplot(2,2,4); imagesc(GroupRGB); title(' Cells outlines on Original Image');
    
end
%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow



% %%% Saves the segmented image, without the splitting individual cells
Var.Img.(Var.Analysis.SegPhiImgOut{CallNum}) = LabelGroupCell;

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;

%%
% %Dilate and erode image to fill cells
% %First Dilation And erosion With large structuring element to maximize the
% %filling of the found cells
% SE = strel('disk',LargeDisk);
% BWDilate_10 = imdilate(BWEdge,SE);
% BWFill = imfill(BWDilate_10,'holes');
% BWErode = imerode(BWFill,SE);
% %resulting BW image is BWLargeArea
% BWLargeArea = BWErode;
% %Label Area found and define properites
% [LabelArea,NumberArea] = bwlabel(BWLargeArea,8);
% props = regionprops(LabelArea,'BoundingBox', 'Area');
% %Second dilation and erosion using a small disk to minimize artifact in the
% %segmentation
% SE = strel('disk',SmallDisk);
% BWDilate_1 = imdilate(BWEdge,SE);
% BWFill = imfill(BWDilate_1,'holes');
% BWErode = imerode(BWFill,SE);
% %Resulting BW image is BWSmallarea
% BWSmallArea = BWErode;
% %Label Area found and define properites
% [LabelSmall,NumberSmall] = bwlabel(BWSmallArea,8);
% Smallprops = regionprops(LabelSmall,'Area', 'PixelIdxList');
%
% %DIlate Edge image with medium structuring element
% SE = strel('disk',MedDisk);
% BWDilate_5 = imdilate(BWEdge,SE);
%
% %Remove from found areas elements which have anb area smaller than the
% %threshold area
% if NumberSmall > 0
%     for n = 1:NumberSmall
%         if Smallprops(n).Area < ThreshArea
%             BWSmallArea(Smallprops(n).PixelIdxList) = 0;
%         end
%     end
% end
% if Debug
%     figure(103);
%     subplot(2,2,1); imagesc(BWSmallArea); title('Small area');
%     subplot(2,2,2); imagesc(BWDilate_5); title('Medium DIlate area');
%     subplot(2,2,3); imagesc(BWLargeArea); title('Large area');
%     subplot(2,2,4); imagesc(LabelArea); title('Label Large area');
% end
% %Define some blank images which will be used further in the analysis
% BWHull = zeros(size(BWErode));
% BWCellSplit = zeros(size(BWErode));
% CellSplit = zeros(size(OrigImage));
% BWHull = zeros(size(OrigImage));
% %Border parameter sets the number of pixels added on each sides of the
% %bounding box defined around the large areas
% Border = 2;
% % takes sequentially every Large area found and tries to split it into
% % smaller cells
% if NumberArea > 0
%     for n = 1:NumberArea
%         if props(n).Area > ThreshArea
%             % Create image containing only one large area
%             %use the bounding box to define the corners of the small image
%             ULCorner = round(props(n).BoundingBox(1:2));
%             Width = props(n).BoundingBox(3:4);
%             Corners = [ULCorner(2),ULCorner(2)+Width(2),ULCorner(1),ULCorner(1)+Width(1) ];
%             %Check if corners are within the size of the full image
%             if ULCorner(1) + Width(1) > size(BWErode,2)
%                 Corners(4) = size(BWErode,2);
%             end
%             if ULCorner(2) + Width(2) > size(BWErode,1)
%                 Corners(2) = size(BWErode,1);
%             end
%             %Cut full images to keep only the small region that contain the
%             %area of interest
%             LabelCut = CutWithBorder(LabelArea, Corners, Border);
%             EdgeCut = CutWithBorder(BWEdge, Corners, Border);
%             SmallCut = CutWithBorder(BWSmallArea, Corners, Border);
%             Dilate_5_Cut = CutWithBorder(BWDilate_5, Corners, Border);
%             % remove other object in the cut image
%             Outside = find(LabelCut ~= n);
%             EdgeCut(Outside) = 0;
%             SmallCut(Outside) = 0;
%             % Set all small object == 1
%             InSmall = find(SmallCut~=0);
%             SmallCut(InSmall) = 1;
%             %Set all edges outside of the small objects to 0
%             OutSmall = find(SmallCut==0);
%             SmallEdgeCut = EdgeCut;
%             SmallEdgeCut(OutSmall) = 0;
%
%             %XY coordinates of all edges found
%             [YEdge, XEdge] = find(EdgeCut == 1);
%             %XY coordinates of edges within small objects
%             [SYEdge, SXEdge] = find(SmallEdgeCut == 1);
%             %Find pixels at the border of a convex shape containing all the
%             %XY coordinates of the edges contained in the small objects
%             if ~isempty(SXEdge)
%             if std(SXEdge) == 0
%                 SXEdge(1) = SXEdge(1) +1;
%             end
%             if std(SYEdge) == 0
%                 SYEdge(1) = SYEdge(1) +1;
%             end
%             end
%             BorderPix = convhull(SXEdge,SYEdge);
%
%             %Find skeleton image of the medium objects
%             SEC_morph = bwmorph(Dilate_5_Cut,'skel',Inf);
%             % Fill Skeleton
%             SEC_fill = imfill(SEC_morph, 'holes');
%             %Subtract filled image by skeleton: Only filled area remain
%             SEC_cell = SEC_fill - SEC_morph;
%             %Calculate distance from filled area
%             SEC_Dist = bwdist(SEC_cell);
%
%             if Debug
%                 figure(108)
%                 subplot(2,2,1); imagesc(SEC_morph); title('SEC_morph');
%                 subplot(2,2,2); imagesc(SEC_fill); title('SEC_fill');
%                 subplot(2,2,3); imagesc(SEC_cell); title('SEC_cell');
%                 subplot(2,2,4); imagesc(SEC_Dist); title('SEC_Dist');
%
%                 AddCut = LabelCut + SmallCut;
%                 AddCut = AddCut + 2.*SEC_morph;
%                 figure(109)
%                 imagesc(AddCut)
%                 hold on
%                 plot(XEdge, YEdge, 'sb')
%                 plot(SXEdge, SYEdge, 'oc')
%                 plot(SXEdge(BorderPix), SYEdge(BorderPix), '-r')
%                 hold off
%
%                 pause(3)
%             end
%             %Further analysis if border pixels where found
%             if ~isempty(BorderPix)
%                 %XY coordinates of Large area
%                 [YIn, XIn] = find(LabelCut == n);
%                 %Find which coordinates of the large area which are
%                 %contained within the convex border of the small area
%                 TestConvHull = inpolygon(XIn,YIn,SXEdge(BorderPix), SYEdge(BorderPix));
%
%                 %Find the indices of the XY coordinates which are within
%                 %the boders
%                 InConvHull = find(TestConvHull == 1);
%                 IndexIN = sub2ind(size(LabelCut),YIn(InConvHull),XIn(InConvHull));
%                 %Define blank images of ther size of the cut images
%                 BWConvHull = zeros(size(LabelCut));
%                 CellSplitCut = zeros(size(LabelCut));
%                 %Set pixels of BWConvHull within the borders to 1
%                 BWConvHull(IndexIN) = 1;
%                 %Split individual cells
%                 %Set SEC_Dist outside of the borders to -inf;
%                 SEC_Dist(BWConvHull==0) = -Inf;
%                 %Define watershed region in SEC Distance image
%                 WatershedCut = watershed(SEC_Dist,4);
%                 %Clear objects touching the border
%                 WatershedCut = imclearborder(WatershedCut,4);
%                 if Debug
%                     figure(111); subplot(2,1,1); imagesc(WatershedCut)
%                 	subplot(2,1,2); imagesc(WatershedCut)
%                 end
%                 %Set Waterched image to black and white
%                 CellSplitCut = im2bw(WatershedCut, 0.5);
%                 %Label found watershed regions
%                 [CellSplitCutLabel, CellNum] = bwlabel(CellSplitCut, 4);
%                 %IF more than one distinct region found
%                 if CellNum>1
%                     Cellprops = regionprops(CellSplitCutLabel,'Area', 'PixelIdxList');
%                     %Remove regions with area smaller than threshold
%                     for m = 1:CellNum
%                         if Cellprops(m).Area < ThreshArea
%                             CellSplitCut(Cellprops(m).PixelIdxList) = 0;
%                         end
%                         %Fill regions if they are contained within another
%                         %region
%                         CellSplitCut = imfill(CellSplitCut,8, 'holes');
%                     end
%                     %Transfer found cells in large image
%                     CellSplit(Corners(1):Corners(2),Corners(3):Corners(4)) = CellSplit(Corners(1):Corners(2),Corners(3):Corners(4)) + CellSplitCut(Border+1:end-Border,Border+1:end-Border);
%                 else
%                     %Transfer found Split cells in large image
%                     CellSplit(Corners(1):Corners(2),Corners(3):Corners(4)) = CellSplit(Corners(1):Corners(2),Corners(3):Corners(4)) +  BWConvHull(Border+1:end-Border,Border+1:end-Border);
%                 end
%
%                 %Transfer found grouped cells in large image
%                 BWHull(Corners(1):Corners(2),Corners(3):Corners(4)) = BWHull(Corners(1):Corners(2),Corners(3):Corners(4)) + BWConvHull(Border+1:end-Border,Border+1:end-Border);
%             end
%         end
%     end
% end
%
% if Debug
%     figure(104); imagesc(CellSplit); title('Split Cells');
%     figure(105); imagesc(BWHull); title('Group Cells');
%      pause(1)
% end
%
% %Test found cells to make sure they are within the specified dimension
% [LabelHull,NumHull] = bwlabel(BWHull);
% Hullprops = regionprops(LabelHull,'Area', 'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList');
% for m = 1:NumHull
%     if Hullprops(m).MajorAxisLength < MinDiam |Hullprops(m).MinorAxisLength > MaxDiam |Hullprops(m).Area < ThreshArea
%         BWHull(Hullprops(m).PixelIdxList) = 0;
%     end
% end
% [LabelHull,NumHull] = bwlabel(BWHull);
%
% [LabelCell,NumCell] = bwlabel(CellSplit);
% Cellprops = regionprops(LabelCell,'Area','MajorAxisLength', 'MinorAxisLength', 'PixelIdxList');
% for m = 1:NumCell
%     if Cellprops(m).MajorAxisLength < MinDiam |Cellprops(m).MinorAxisLength > MaxDiam | Cellprops(m).Area < ThreshArea
%         CellSplit(Cellprops(m).PixelIdxList) = 0;
%     end
% end
% [LabelCell,NumCell] = bwlabel(CellSplit);
%
% if Debug
%     figure(107); imagesc(LabelHull); title('Labelled Hull')
%     figure(106); imagesc(LabelCell); title('Labelled Cell')
% end
%
% %Remove objects touching border if necessary
% if strncmpi(IncludeEdge,'N',1) == 1
%     LabelHull = imclearborder(LabelHull,8);
%     LabelCell = imclearborder(LabelCell,8);
% end
% [LabelHull,NumHull] = bwlabel(BWHull);
% [LabelCell,NumCell] = bwlabel(CellSplit);
% FinalLabelMatrixImage = LabelCell;
%
%
% %%%%%%%%%%%%%%%%%%%%%%
% %%% DISPLAY RESULTS %%%
% %%%%%%%%%%%%%%%%%%%%%%
% drawnow
%
%
%
% if strcmp(Var.Figure.Display, 'on')
%     if sum(sum(FinalLabelMatrixImage)) >= 1
%         ColoredLabelMatrixImage = label2rgb(FinalLabelMatrixImage);
%     else
%         ColoredLabelMatrixImage = FinalLabelMatrixImage;
%     end
%
%     %%% Calculates the object outlines, which are overlaid on the original
%     %%% image and displayed in figure subplot (2,2,4).
%     %%% Creates the structuring element that will be used for dilation.
%     StructuringElement = strel('square',3);
%     %%% Converts the FinalLabelMatrixImage to binary.
%     FinalBinaryImage = im2bw(FinalLabelMatrixImage,0.5);
%     %%% Dilates the FinalBinaryImage by one pixel (8 neighborhood).
%     DilatedBinaryImage = imdilate(FinalBinaryImage, StructuringElement);
%     %%% Subtracts the FinalBinaryImage from the DilatedBinaryImage,
%     %%% which leaves the PrimaryObjectOutlines.
%     PrimaryObjectOutlines = DilatedBinaryImage - FinalBinaryImage;
%     %%% Overlays the object outlines on the original image.
%     ObjectOutlinesOnOrigImage = OrigImage;
%     %%% Determines the grayscale intensity to use for the cell outlines.
%     LineIntensity = max(OrigImage(:));
%     ObjectOutlinesOnOrigImage(PrimaryObjectOutlines == 1) = LineIntensity;
%
%     figure(find((strcmp(Var.Figure.List, 'SegYP'))))
%     %%% A subplot of the figure window is set to display the original image.
%     subplot(2,2,1); imagesc(OrigImage); title('Input Image');
%     %%% A subplot of the figure window is set to display the colored label
%     %%% matrix image.
%     subplot(2,2,2); imagesc(ColoredLabelMatrixImage); title('Segmented Cells');
%     %%% A subplot of the figure window is set to display the Overlaid image,
%     %%% where the maxima are imposed on the inverted original image
%     subplot(2,2,3); imagesc(EnhancedInvertedImage); title('Inverted enhanced contrast image');
%     %%% A subplot of the figure window is set to display the inverted original
%     %%% image with watershed lines drawn to divide up clusters of objects.
%     subplot(2,2,4); imagesc(ObjectOutlinesOnOrigImage); title('Cell Outlines on Input Image');
%
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% SAVE DATA TO HANDLES STRUCTURE %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drawnow
%
%
%
% % %%% Saves the segmented image, without the splitting individual cells
% % Var.Img.GroupSegYeast = LabelHull;
% %
% % %%% Saves the final segmented label matrix image to the handles structure.
% % Var.Img.SplitSegYeast = FinalLabelMatrixImage;
% Var.Img.SplitCell = LabelHull;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Sub-function          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ImOut = CutWithBorder(ImIn, Corners, Border);
%corner: [Ymin, Ymax, Xmin, Xmax]
ImMid = ImIn(Corners(1):Corners(2),Corners(3):Corners(4));

ImOut = zeros(size(ImMid,1)+2*Border,size(ImMid,2)+2*Border);

ImOut(Border+1:size(ImMid,1)+Border,Border+1:size(ImMid,2)+Border)  = ImMid;







