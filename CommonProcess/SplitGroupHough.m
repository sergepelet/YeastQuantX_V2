function Var = SplitGroupHough(Var, CallNum)
tic
Debug = 0;      %Set to 1 to display segmentation images and 0 not to
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
MaxDiam = Var.SegmentationParameters(ParaNum).MaxDiameter;



%%
%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%

% Get images from structure
SplitName = Var.Analysis.SplitCellOut{CallNum};
GroupName = Var.Analysis.SplitGroupIn{CallNum};
GroupImg = Var.Img.(GroupName);
OrigImg = Var.Img.(Var.Analysis.SplitGroupPhiImg{CallNum});
OrigImg =  OrigImg - mean(OrigImg(:));

% if second bright field image is available use it
if isfield(Var.Analysis, 'SplitGroupPhiSecImg') && ~isempty(Var.Analysis.SplitGroupPhiSecImg)
    SecImg = double(Var.Img.(Var.Analysis.SegPhiSecImg{CallNum}));
    SecZ = 1;
    SecImg =  SecImg - mean(SecImg(:));
    %Add first and second BF image
    %OrigImg = SecImg + OrigImg;
else SecZ = 0;
end

if Debug
    figure(102)
    subplot(2,2,1);imagesc(OrigImg); title('OrigImg')
    subplot(2,2,2);imagesc(GroupImg); title('GroupImg')
    if SecZ
        subplot(2,2,3);imagesc(SecImg); title('SecImg')
    end
end

%% Analysie all Groups sequentially
%Label Group image
[LabelGroupCell,NumGroupCell] = bwlabel(GroupImg);
Groupprops = regionprops(LabelGroupCell,'BoundingBox','ConvexArea','ConvexImage', 'Area');

% %Get Label from tracking
% FrameIter = Var.Analysis.FrameIter;
% if isfield(Var, 'Measurements')
%     Label = Var.Measurements.(GroupName).Label{FrameIter};
% else
    Label = [1:NumGroupCell];
%end

%Border parameter sets the number of pixels added on each sides of the
%bounding box defined around the large areas
Border = 10;
FinalSplitCell  = zeros(size(GroupImg));
PrevAddedSplitCell  = zeros(size(GroupImg));
TotNumSplitCells = 0;
if NumGroupCell > 0
    for n = 1:NumGroupCell
        % Create image containing only one large area
        %use the bounding box to define the corners of the small image
        ULCorner = round(Groupprops(n).BoundingBox(1:2));
        Width = Groupprops(n).BoundingBox(3:4);
        if Width(1) < 32-2*Border
            Width(1) = 32-2*Border;
        end
        if Width(2) < 32-2*Border
            Width(2) = 32-2*Border;
        end
        Corners = [ULCorner(2),ULCorner(2)+Width(2),ULCorner(1),ULCorner(1)+Width(1) ];
        %Check if corners are within the size of the full image
        if ULCorner(1) + Width(1) > size(GroupImg,2)
            Corners(4) = size(GroupImg,2);
        end
        if ULCorner(2) + Width(2) > size(GroupImg,1)
            Corners(2) = size(GroupImg,1);
        end

        %Cut full images to keep only the small region that contain the
        %area of interest
        CutGroup = CutWithBorder(LabelGroupCell, Corners, Border);
        CutOrig = CutWithBorder(OrigImg, Corners, Border);
        

        %Remove Pixels belonging to other objects
        CutGroup(CutGroup~= n) = 0;
        CutOrig(CutGroup~= n) = 0;
        if min(size(CutOrig)) < 32
            fprintf('Small Cut')
            %Skip to the analysis of the next group of cells
            continue
        end
        %Calculate Hough transform to find round objects in image
        [AccumImg, CellCenter, CellRad] = CircularHough_Grd(CutOrig, [MinDiam/2,MaxDiam/2],  10, 20, 1); %15, 25

        if SecZ
            CutSec = CutWithBorder(SecImg, Corners, Border);
            CutSec(CutGroup~= n) = 0;
            [SecAccumImg, SecCellCenter, SecCellRad] = CircularHough_Grd(CutSec, [MinDiam/2,MaxDiam/2],  10, 20, 1); %15, 25
            CellCenter = [CellCenter;SecCellCenter];
            CellRad = [CellRad;SecCellRad];
        end



        %Get pixels within objject
        [YGr,XGr] = find(CutGroup == n);
        %Loop through all cells in group
        NbCell =length(CellRad);
        ImgCell = zeros(size(CutGroup));
        CenterCell = zeros(size(CutGroup));
        
        %Sort Radius to start with larger cells
        [CellRad, RadOrder] = sort(CellRad, 'descend');
        CellCenter = CellCenter(RadOrder,:);
        
        IterCell = 0;
        for C = 1:NbCell

            % find pixels within circle given by radius and center
            Angle = [0 : (2 * pi / 100) : (2 * pi)];
            XPerim = CellRad(C) * cos(Angle) + CellCenter(C,1);
            YPerim = CellRad(C) * sin(Angle) + CellCenter(C,2);
            InCell = inpolygon(XGr,YGr,XPerim,YPerim);
            %Find indices within circle
            InInd = sub2ind(size(CutGroup),YGr(InCell),XGr(InCell));
            %Compare circel and object area
            CircleArea = pi*CellRad(C)^2;
            OverlapArea = length(InInd);
            %If less than 80% of full circle, than ignore this object
            %Meaning there is low overlap of hough circle with segmented
            %group of cells
            if OverlapArea > CircleArea*0.8
                %Check if the center is not in a previously identified Circle 
                if ImgCell(round(CellCenter(C,2)), round(CellCenter(C,1))) == 0
                    %Count number of cells
                    IterCell = IterCell+1;
                    %Set pixels of ImgCell to IterCell
                    ImgCell(InInd) = IterCell;
                    CenterCell(round(CellCenter(C,2)), round(CellCenter(C,1))) = 1;
                else
                    %If previously identified find cell number
                    SameCell = ImgCell(round(CellCenter(C,2)), round(CellCenter(C,1)));
                    %And Add pixels from this new circle to that given cell 
                    ImgCell(InInd) = SameCell;
                end
            end
        end


        %Label ImgCell to find cells which overlapp
        [LabelImgCell,NumImgCell] = bwlabel(ImgCell);
        %Loop through all lebeled cells
        for C = 1:NumImgCell
            %Find Number of cells contained in each region
            NumCircles = find(unique(ImgCell(LabelImgCell == C)));
            NumCircles = NumCircles(NumCircles > 0);
            %if more than one region Split by watershed
            if length(NumCircles) > 1
                %Calculate dist transform based on centers
                DistCenter = CenterCell;
                DistCenter(LabelImgCell~= C) = 0;
                DistCenter = bwdist(DistCenter);
                %Calculate distance transform based on cell shape
                DistCell = zeros(size(DistCenter));
                DistCell(LabelImgCell ~= C) = 1;
                DistCell = -1.*bwdist(DistCell);
                %Add 2 dist transform to get watershed regions
                SumDist = DistCell +DistCenter; %./2;
                %SumDist(SumDist<-15) = -15;
                %Set Pixels outside cells to -inf
                SumDist(LabelImgCell ~= C) = -Inf;
                %Supresse reginonal minima
                SumDistRegMin = imhmin(SumDist,3);
                %Calculate watershed regions
                WaterCell = double(watershed(SumDistRegMin));
                WaterCell(LabelImgCell ~= C) = -1;
                %Set watershed borders to 0 in ImgCell image
                ImgCell(WaterCell==0) = 0;
                if Debug
                    figure(101)
                    subplot(2,2,1); imagesc(SumDist); title('SumDist')
                    RegionalMin = imregionalmin(SumDistRegMin);
                    subplot(2,2,2); imagesc(RegionalMin); title('RegionalMin')
                    subplot(2,2,3); imagesc(WaterCell); title('WaterCell')

                    subplot(2,2,4); imagesc(ImgCell); title('ImgCell')
                    pause(0)
                end
            end
        end
        %Label again ImgCell
        [LabelImgCell,NumImgCell] = bwlabel(ImgCell);
        

        %Check if more cell found 
        if NumImgCell > IterCell 
            %Get area from Cells
            CellProps = regionprops(LabelImgCell,'Area','PixelIdxList');
            %Sort areas
            [AllArea, AreaIndex] =  sort(cell2mat({CellProps(:).Area}));
            %For all the smallest areas until the cell number matches set
            %to zero
            for i = 1:(NumImgCell-IterCell)
                ImgCell(CellProps(AreaIndex(i)).PixelIdxList) = 0;
            end
        elseif NumImgCell < IterCell 
                %Supresse reginonal minima
                SumDistRegMin = imhmin(SumDist,1);
                %Calculate watershed regions
                WaterCell = watershed(SumDistRegMin);
                WaterCell(LabelImgCell~= C) = -1;
                %Set watershed borders to 0 in ImgCell image
                ImgCell(WaterCell==0) = 0;
        end

        %Label again image
        [LabelImgCell,NumImgCell] = bwlabel(ImgCell);
        %Check again if cell count matche number of cells found
        if IterCell ~= NumImgCell
            fprintf(['Still Pb With Cell Splitting:, ',num2str(IterCell), 'Found Cells and ' num2str(NumImgCell), 'Labelled Cells'])
            warning(['PB with cell Split in SplitGroupHough in analysis of ', Var.Analysis.CurrentFolder,' at frame No ', num2str(Var.Analysis.CurrentFrame)] )
            FileName = [Var.Analysis.OutPath,'_F', num2str(Var.Analysis.CurrentFrame),'_SGH1.mat'];
            save(FileName, 'Var');
        end

        
        %Transfer group cell in Final image
        FinalSplitCell(Corners(1):Corners(2),Corners(3):Corners(4)) = FinalSplitCell(Corners(1):Corners(2),Corners(3):Corners(4)) ...
            + LabelImgCell(Border+1:end-Border,Border+1:end-Border);

        %% Check Number of cells in group compared to previous frame
%             if strcmp(Var.Analysis.FrameAnalyzed, 'First') || strcmp(Var.Analysis.FrameAnalyzed, 'Single')
%                 %For first frame, Write cell number to measurement structure
%                 Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%             else
%                 %Get number of cells in group from Measurment VAR
%                 PrevLabelIndex = find(Var.Measurements.(GroupName).Label{FrameIter-1} == Label(n));
% 
%                 if isempty(PrevLabelIndex) || Label(n) > length(Var.Measurements.(GroupName).NumCell{FrameIter-1})
%                     %If there is no corresponding label in previous image only
%                     %assign the number of cells in var
%                     Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%                 else
%                     %If group of cell was already present, it has a label
%                     %Assign number of cells in group to a previousNumCell
%                     PreviousNumCell = Var.Measurements.(GroupName).NumCell{FrameIter-1}(Label(n));
%                     %get area by counting number of pixels in area
% %                     PrevArea = length(Var.Measurements.(GroupName).PixelList{FrameIter-1,PrevLabelIndex});
% %                     NewArea = length(Var.Measurements.(GroupName).PixelList{FrameIter,n});
%                     %Check if we have more cells now
%                     if PreviousNumCell < NumImgCell
%                         Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%                         %NOT USED FOR NOW Hough transform should not
%                         %over segment
%                         %Control if area has increased
%                         %If the previous area + threshold area is smaller
%                         %than the new area, One can think that a new cell is
%                         %present
%                         %                           fprintf(['GroupNum: ', num2str(n), ' OldCellNum: ',num2str(PreviousNumCell), ' NewCellNum: ', num2str(NumImgCell),' OldArea: ', num2str(PrevArea), ' NewArea: ', num2str(NewArea)])
%                         %                             if PrevArea + ThreshArea < NewArea
%                         %                                 Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%                         %                                 RedoCellSplit = 0;
%                         %                                 fprintf(' No Split\n')
%                         %                             else
%                         %                                 RedoCellSplit = 1;
%                         %                                 fprintf(' Split\n')
%                         %                             end
%                         %Check if we have less cells now
%                     elseif PreviousNumCell > NumImgCell
%                         % fprintf(['GroupNum: ', num2str(n), ' OldCellNum: ',num2str(PreviousNumCell), ' NewCellNum: ', num2str(NumImgCell),' OldArea: ', num2str(PrevArea), ' NewArea: ', num2str(NewArea)])
%                         %Cells should not diseapear but if there is flow it can be
%                         %possible. Check overlap between Cells
%                         PrevCellLabels = find(Var.Measurements.(SplitName).([GroupName, 'Label']){FrameIter-1} == (Label(n)));
%                         %Make image with cells present previously in the
%                         %group
%                         PrevCellImg = zeros(size(GroupImg));
%                         for PC = 1:length(PrevCellLabels)
%                             PrevCellImg(Var.Measurements.(SplitName).PixelList{FrameIter-1, PrevCellLabels(PC)}) = PC;
%                         end
%                         %Cut image to make it smaller
%                         CompareCells = CutWithBorder(PrevCellImg, Corners, Border);
% 
%                         %Set Pixels recognized by new segmentation to zero
%                         CompareCells(LabelImgCell > 0) = 0;
%                         %Set Pixels outside of current group to 0
%                         CompareCells(CutGroup == 0) = 0;
% 
%                         %                             [LabelCompareCells,NumCompCell] = bwlabel(CompareCells);
%                         %                             CompProps = regionprops(LabelCompareCells,'PixelIdxList', 'Area');
% 
%                         NumAddedCell = 0;
%                         for PC = 1:length(PrevCellLabels)
%                             %Check number of pixels remaining from previous
%                             %cell
%                             NoOverlapPix = find(CompareCells == PC);
%                             SizeOverlap = length(NoOverlapPix);
%                             %Get size from current group also present in
%                             %the no overlap pixels
%                             SizeFromGroup = length(find(CutGroup(NoOverlapPix) == 1));
%                             %Get from the previous cell
%                             SizePrevCell = length(Var.Measurements.(SplitName).PixelList{FrameIter-1, PrevCellLabels(PC)});
%                             %Compare overlap with previous cell size and
%                             %total size from group
%                             if SizeOverlap>0.8*SizePrevCell && SizeOverlap > 0.5*SizeFromGroup
%                                 PrevAddedSplitCell(Var.Measurements.(SplitName).PixelList{FrameIter-1, PrevCellLabels(PC)}) = max(PrevAddedSplitCell(:))+1;
%                                 NumAddedCell = NumAddedCell+1;
%                             end
% 
%                         end
%                         NumImgCell = NumImgCell + NumAddedCell;
%                         Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%                         if Debug > 0
%                             figure(103)
%                             subplot(2,2,1);imagesc(CompareCells); title('CutPrevCell')
%                             subplot(2,2,2);imagesc(LabelImgCell); title('LabelImgCell')
%                             %DiffImg = LabelImgCell - CutPrevCell;
%                             subplot(2,2,3);imagesc(CutGroup); title('CutGroup')
%                             subplot(2,2,4);imagesc(PrevAddedSplitCell); title('PrevAddedSplitCell')
%                         end
%                     elseif PreviousNumCell == NumImgCell
%                         Var.Measurements.(GroupName).NumCell{FrameIter}(Label(n)) = NumImgCell;
%                     end
%                 end
%             end
            % Calculate total number of cells in image
            TotNumSplitCells = TotNumSplitCells + NumImgCell;
    end

    AllCells = FinalSplitCell + PrevAddedSplitCell;
    [LabelSplitCell,NumGroupCell] = bwlabel(AllCells);
    if NumGroupCell ~= TotNumSplitCells
        OrigFinalSplitCell = FinalSplitCell;
        OrigFinalSplitCell(OrigFinalSplitCell>0) =1;
        SE = strel('disk',2);
        OrigFinalSplitCell = imdilate(OrigFinalSplitCell,SE);
        AllAddedCell = max(PrevAddedSplitCell(:));
        PrevAddedSplitCell(OrigFinalSplitCell>0) = 0;
        [PrevAddedSplitCell,NumAddedCell] = bwlabel(PrevAddedSplitCell);
        if NumAddedCell > AllAddedCell
            AddedProps = regionprops(PrevAddedSplitCell,'PixelIdxList', 'Area');
            [AllArea, AreaIndex] =  sort(cell2mat({AddedProps(:).Area}));
            for i = 1:(NumAddedCell-AllAddedCell)
                PrevAddedSplitCell(AddedProps(AreaIndex(i)).PixelIdxList) = 0;
            end
        end

        FinalSplitCell = FinalSplitCell + PrevAddedSplitCell;
        [LabelSplitCell,NumGroupCell] = bwlabel(FinalSplitCell);
        if NumGroupCell ~= TotNumSplitCells
            if Debug > 0
                figure(104)
                subplot(2,2,1);imagesc(FinalSplitCell); title('FinalSplitCell')
                subplot(2,2,2);imagesc(PrevAddedSplitCell); title('PrevAddedSplitCell')
                subplot(2,2,4);imagesc(OrigFinalSplitCell); title('OrigFinalSplitCell')
            end
            warning(['Added cell overlap?? in analysis of ', Var.Analysis.CurrentFolder,' at frame No ', num2str(Var.Analysis.CurrentFrame)] )
           	FileName = [Var.Analysis.OutPath,'_F', num2str(Var.Analysis.CurrentFrame),'_SGH2.mat'];
            save(FileName, 'Var');
        end
    end
end

%Label final image
[LabelSplitCell,NumGroupCell] = bwlabel(FinalSplitCell);
%%

%Create Diameter map
tmp = regionprops(LabelSplitCell,'EquivDiameter');
Diameters = [0;cat(1,tmp.EquivDiameter)];

%%% Locate objects with diameter smalle than min diameter
    %%% Create image with object intensity equal to the diameter
    DiameterMap = Diameters(LabelSplitCell+1);
    %%% Remove objects that are too small
    LabelSplitCell(DiameterMap < MinDiam*.5) = 0;
%ReLabel final image
[LabelSplitCell,NumGroupCell] = bwlabel(LabelSplitCell);

%%

%Display
if strcmp(Var.Figure.Display, 'on')

    FigNum = find(strcmp(Var.Figure.List, 'SplitCell'));
    figure(FigNum(CallNum))

    GroupR = (OrigImg-min(OrigImg(:)))./(max(OrigImg(:))-min(OrigImg(:)));
    GroupR(bwperim(LabelSplitCell)) = 1;
    GroupG = (OrigImg-min(OrigImg(:)))./(max(OrigImg(:))-min(OrigImg(:)));
    %GroupG(BWLowGroupCell) = 0.5;
    GroupRGB(:,:,1) = GroupR;
    GroupRGB(:,:,2) = GroupG;
    GroupRGB(:,:,3) = (OrigImg-min(OrigImg(:)))./(max(OrigImg(:))-min(OrigImg(:)));
    subplot(2,2,1); imagesc(GroupImg), title('GroupImg')
    subplot(2,2,2); imagesc(OrigImg), title('OrigImg')
    subplot(2,2,3); imagesc(LabelSplitCell), title(['LabelSplitCell', num2str(NumGroupCell), 'Num Cells'] )

    subplot(2,2,4); image(GroupRGB), title('OrigImg With borders')


end

%%

%Save to Var
Var.Img.(Var.Analysis.SplitCellOut{CallNum}) = LabelSplitCell;
%Save X-Y center position of objects
ObjProps = regionprops(LabelSplitCell,'Centroid');
if length(ObjProps) > 0
    for i = 1:length(ObjProps)
        Var.Measurements.(Var.Analysis.SplitCellOut{CallNum}).CenterX(i,1) = ObjProps(i).Centroid(1);
        Var.Measurements.(Var.Analysis.SplitCellOut{CallNum}).CenterY(i,1) = ObjProps(i).Centroid(2);
        Var.Measurements.(Var.Analysis.SplitCellOut{CallNum}).SegLabel(i,1) = i;
    end
else
    error( 'No object found in in Image by Segmentation')
end
%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Sub-function          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [accum, varargout] = CircularHough_Grd(img, radrange, varargin)
%Detect circular shapes in a grayscale image. Resolve their center
%positions and radii.
%
%  [accum, circen, cirrad, dbg_LMmask] = CircularHough_Grd(
%      img, radrange, grdthres, fltr4LM_R, multirad, fltr4accum)
%  Circular Hough transform based on the gradient field of an image.
%  NOTE:    Operates on grayscale images, NOT B/W bitmaps.
%           NO loops in the implementation of Circular Hough transform,
%               which means faster operation but at the same time larger
%               memory consumption.
%
%%%%%%%% INPUT: (img, radrange, grdthres, fltr4LM_R, multirad, fltr4accum)
%
%  img:         A 2-D grayscale image (NO B/W bitmap)
%
%  radrange:    The possible minimum and maximum radii of the circles
%               to be searched, in the format of
%               [minimum_radius , maximum_radius]  (unit: pixels)
%               **NOTE**:  A smaller range saves computational time and
%               memory.
%
%  grdthres:    (Optional, default is 10, must be non-negative)
%               The algorithm is based on the gradient field of the
%               input image. A thresholding on the gradient magnitude
%               is performed before the voting process of the Circular
%               Hough transform to remove the 'uniform intensity'
%               (sort-of) image background from the voting process.
%               In other words, pixels with gradient magnitudes smaller
%               than 'grdthres' are NOT considered in the computation.
%               **NOTE**:  The default parameter value is chosen for
%               images with a maximum intensity close to 255. For cases
%               with dramatically different maximum intensities, e.g.
%               10-bit bitmaps in stead of the assumed 8-bit, the default
%               value can NOT be used. A value of 4% to 10% of the maximum
%               intensity may work for general cases.
%
%  fltr4LM_R:   (Optional, default is 8, minimum is 3)
%               The radius of the filter used in the search of local
%               maxima in the accumulation array. To detect circles whose
%               shapes are less perfect, the radius of the filter needs
%               to be set larger.
%
% multirad:     (Optional, default is 0.5)
%               In case of concentric circles, multiple radii may be
%               detected corresponding to a single center position. This
%               argument sets the tolerance of picking up the likely
%               radii values. It ranges from 0.1 to 1, where 0.1
%               corresponds to the largest tolerance, meaning more radii
%               values will be detected, and 1 corresponds to the smallest
%               tolerance, in which case only the "principal" radius will
%               be picked up.
%
%  fltr4accum:  (Optional. A default filter will be used if not given)
%               Filter used to smooth the accumulation array. Depending
%               on the image and the parameter settings, the accumulation
%               array built has different noise level and noise pattern
%               (e.g. noise frequencies). The filter should be set to an
%               appropriately size such that it's able to suppress the
%               dominant noise frequency.
%
%%%%%%%% OUTPUT: [accum, circen, cirrad, dbg_LMmask]
%
%  accum:       The result accumulation array from the Circular Hough
%               transform. The accumulation array has the same dimension
%               as the input image.
%
%  circen:      (Optional)
%               Center positions of the circles detected. Is a N-by-2
%               matrix with each row contains the (x, y) positions
%               of a circle. For concentric circles (with the same center
%               position), say k of them, the same center position will
%               appear k times in the matrix.
%
%  cirrad:      (Optional)
%               Estimated radii of the circles detected. Is a N-by-1
%               column vector with a one-to-one correspondance to the
%               output 'circen'. A value 0 for the radius indicates a
%               failed detection of the circle's radius.
%
%  dbg_LMmask:  (Optional, for debugging purpose)
%               Mask from the search of local maxima in the accumulation
%               array.
%
%%%%%%%%% EXAMPLE #0:
%  rawimg = imread('TestImg_CHT_a2.bmp');
%  tic;
%  [accum, circen, cirrad] = CircularHough_Grd(rawimg, [15 60]);
%  toc;
%  figure(1); imagesc(accum); axis image;
%  title('Accumulation Array from Circular Hough Transform');
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  hold on;
%  plot(circen(:,1), circen(:,2), 'r+');
%  for k = 1 : size(circen, 1),
%      DrawCircle(circen(k,1), circen(k,2), cirrad(k), 32, 'b-');
%  end
%  hold off;
%  title(['Raw Image with Circles Detected ', ...
%      '(center positions and radii marked)']);
%  figure(3); surf(accum, 'EdgeColor', 'none'); axis ij;
%  title('3-D View of the Accumulation Array');
%
%  COMMENTS ON EXAMPLE #0:
%  Kind of an easy case to handle. To detect circles in the image whose
%  radii range from 15 to 60. Default values for arguments 'grdthres',
%  'fltr4LM_R', 'multirad' and 'fltr4accum' are used.
%
%%%%%%%%% EXAMPLE #1:
%  rawimg = imread('TestImg_CHT_a3.bmp');
%  tic;
%  [accum, circen, cirrad] = CircularHough_Grd(rawimg, [15 60], 10, 20);
%  toc;
%  figure(1); imagesc(accum); axis image;
%  title('Accumulation Array from Circular Hough Transform');
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  hold on;
%  plot(circen(:,1), circen(:,2), 'r+');
%  for k = 1 : size(circen, 1),
%      DrawCircle(circen(k,1), circen(k,2), cirrad(k), 32, 'b-');
%  end
%  hold off;
%  title(['Raw Image with Circles Detected ', ...
%      '(center positions and radii marked)']);
%  figure(3); surf(accum, 'EdgeColor', 'none'); axis ij;
%  title('3-D View of the Accumulation Array');
%
%  COMMENTS ON EXAMPLE #1:
%  The shapes in the raw image are not very good circles. As a result,
%  the profile of the peaks in the accumulation array are kind of
%  'stumpy', which can be seen clearly from the 3-D view of the
%  accumulation array. (As a comparison, please see the sharp peaks in
%  the accumulation array in example #0) To extract the peak positions
%  nicely, a value of 20 (default is 8) is used for argument 'fltr4LM_R',
%  which is the radius of the filter used in the search of peaks.
%
%%%%%%%%% EXAMPLE #2:
%  rawimg = imread('TestImg_CHT_b3.bmp');
%  fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, circen, cirrad] = CircularHough_Grd(imgfltrd, [15 80], 8, 10);
%  toc;
%  figure(1); imagesc(accum); axis image;
%  title('Accumulation Array from Circular Hough Transform');
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  hold on;
%  plot(circen(:,1), circen(:,2), 'r+');
%  for k = 1 : size(circen, 1),
%      DrawCircle(circen(k,1), circen(k,2), cirrad(k), 32, 'b-');
%  end
%  hold off;
%  title(['Raw Image with Circles Detected ', ...
%      '(center positions and radii marked)']);
%
%  COMMENTS ON EXAMPLE #2:
%  The circles in the raw image have small scale irregularities along
%  the edges, which could lead to an accumulation array that is bad for
%  local maxima detection. A 5-by-5 filter is used to smooth out the
%  small scale irregularities. A blurred image is actually good for the
%  algorithm implemented here which is based on the image's gradient
%  field.
%
%%%%%%%%% EXAMPLE #3:
%  rawimg = imread('TestImg_CHT_c3.bmp');
%  fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, circen, cirrad] = ...
%      CircularHough_Grd(imgfltrd, [15 105], 8, 10, 0.7);
%  toc;
%  figure(1); imagesc(accum); axis image;
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  hold on;
%  plot(circen(:,1), circen(:,2), 'r+');
%  for k = 1 : size(circen, 1),
%      DrawCircle(circen(k,1), circen(k,2), cirrad(k), 32, 'b-');
%  end
%  hold off;
%  title(['Raw Image with Circles Detected ', ...
%      '(center positions and radii marked)']);
%
%  COMMENTS ON EXAMPLE #3:
%  Similar to example #2, a filtering before circle detection works for
%  noisy image too. 'multirad' is set to 0.7 to eliminate the false
%  detections of the circles' radii.
%
%%%%%%%%% BUG REPORT:
%  This is a beta version. Please send your bug reports, comments and
%  suggestions to pengtao@glue.umd.edu . Thanks.
%
%
%%%%%%%%% INTERNAL PARAMETERS:
%  The INPUT arguments are just part of the parameters that are used by
%  the circle detection algorithm implemented here. Variables in the code
%  with a prefix 'prm_' in the name are the parameters that control the
%  judging criteria and the behavior of the algorithm. Default values for
%  these parameters can hardly work for all circumstances. Therefore, at
%  occasions, the values of these INTERNAL PARAMETERS (parameters that
%  are NOT exposed as input arguments) need to be fine-tuned to make
%  the circle detection work as expected.
%  The following example shows how changing an internal parameter could
%  influence the detection result.
%  1. Change the value of the internal parameter 'prm_LM_LoBndRa' to 0.4
%     (default is 0.2)
%  2. Run the following matlab code:
%     fltr4accum = [1 2 1; 2 6 2; 1 2 1];
%     fltr4accum = fltr4accum / sum(fltr4accum(:));
%     rawimg = imread('Frame_0_0022_portion.jpg');
%     tic;
%     [accum, circen] = CircularHough_Grd(rawimg, ...
%         [4 14], 10, 4, 0.5, fltr4accum);
%     toc;
%     figure(1); imagesc(accum); axis image;
%     title('Accumulation Array from Circular Hough Transform');
%     figure(2); imagesc(rawimg); colormap('gray'); axis image;
%     hold on; plot(circen(:,1), circen(:,2), 'r+'); hold off;
%     title('Raw Image with Circles Detected (center positions marked)');
%  3. See how different values of the parameter 'prm_LM_LoBndRa' could
%     influence the result.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: Beta        Revision: Mar. 07, 2007


%%%%%%%% Arguments and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validation of arguments
if ndims(img) ~= 2 || ~isnumeric(img),
    error('CircularHough_Grd: ''img'' has to be 2 dimensional');
end
if ~all(size(img) >= 32),
    error('CircularHough_Grd: ''img'' has to be larger than 32-by-32');
end

if numel(radrange) ~= 2 || ~isnumeric(radrange),
    error(['CircularHough_Grd: ''radrange'' has to be ', ...
        'a two-element vector']);
end
prm_r_range = sort(max( [0,0;radrange(1),radrange(2)] ));

% Parameters (default values)
prm_grdthres = 10;
prm_fltrLM_R = 8;
prm_multirad = 0.5;
func_compu_cen = true;
func_compu_radii = true;

% Validation of arguments
vap_grdthres = 1;
if nargin > (1 + vap_grdthres),
    if isnumeric(varargin{vap_grdthres}) && ...
            varargin{vap_grdthres}(1) >= 0,
        prm_grdthres = varargin{vap_grdthres}(1);
    else
        error(['CircularHough_Grd: ''grdthres'' has to be ', ...
            'a non-negative number']);
    end
end

vap_fltr4LM = 2;    % filter for the search of local maxima
if nargin > (1 + vap_fltr4LM),
    if isnumeric(varargin{vap_fltr4LM}) && varargin{vap_fltr4LM}(1) >= 3,
        prm_fltrLM_R = varargin{vap_fltr4LM}(1);
    else
        error(['CircularHough_Grd: ''fltr4LM_R'' has to be ', ...
            'larger than or equal to 3']);
    end
end

vap_multirad = 3;
if nargin > (1 + vap_multirad),
    if isnumeric(varargin{vap_multirad}) && ...
            varargin{vap_multirad}(1) >= 0.1 && ...
            varargin{vap_multirad}(1) <= 1,
        prm_multirad = varargin{vap_multirad}(1);
    else
        error(['CircularHough_Grd: ''multirad'' has to be ', ...
            'within the range [0.1, 1]']);
    end
end

vap_fltr4accum = 4; % filter for smoothing the accumulation array
if nargin > (1 + vap_fltr4accum),
    if isnumeric(varargin{vap_fltr4accum}) && ...
            ndims(varargin{vap_fltr4accum}) == 2 && ...
            all(size(varargin{vap_fltr4accum}) >= 3),
        fltr4accum = varargin{vap_fltr4accum};
    else
        error(['CircularHough_Grd: ''fltr4accum'' has to be ', ...
            'a 2-D matrix with a minimum size of 3-by-3']);
    end
else
    % Default filter (5-by-5)
    fltr4accum = ones(5,5);
    fltr4accum(2:4,2:4) = 2;
    fltr4accum(3,3) = 6;
end

func_compu_cen = ( nargout > 1 );
func_compu_radii = ( nargout > 2 );

% Reserved parameters
dbg_on = false;      % debug information
dbg_bfigno = 4;
if nargout > 3,  dbg_on = true;  end


%%%%%%%% Building accumulation array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the image to single if it is not of
% class float (single or double)
img_is_double = isa(img, 'double');
if ~(img_is_double || isa(img, 'single')),
    imgf = single(img);
end

% Compute the gradient and the magnitude of gradient
if img_is_double,
    [grdx, grdy] = gradient(img);
else
    [grdx, grdy] = gradient(imgf);
end
grdmag = sqrt(grdx.^2 + grdy.^2);


% Get the linear indices, as well as the subscripts, of the pixels
% whose gradient magnitudes are larger than the given threshold
grdmasklin = find(grdmag > prm_grdthres);
[grdmask_IdxI, grdmask_IdxJ] = ind2sub(size(grdmag), grdmasklin);

% Compute the linear indices (as well as the subscripts) of
% all the votings to the accumulation array.
% The Matlab function 'accumarray' accepts only double variable,
% so all indices are forced into double at this point.
% A row in matrix 'lin2accum_aJ' contains the J indices (into the
% accumulation array) of all the votings that are introduced by a
% same pixel in the image. Similarly with matrix 'lin2accum_aI'.
rr_4linaccum = double( prm_r_range );
linaccum_dr = [ (-rr_4linaccum(2) + 0.5) : -rr_4linaccum(1) , ...
    (rr_4linaccum(1) + 0.5) : rr_4linaccum(2) ];

lin2accum_aJ = floor( ...
    double(grdx(grdmasklin)./grdmag(grdmasklin)) * linaccum_dr + ...
    repmat( double(grdmask_IdxJ)+0.5 , [1,length(linaccum_dr)] ) ...
    );
lin2accum_aI = floor( ...
    double(grdy(grdmasklin)./grdmag(grdmasklin)) * linaccum_dr + ...
    repmat( double(grdmask_IdxI)+0.5 , [1,length(linaccum_dr)] ) ...
    );

% Clip the votings that are out of the accumulation array
mask_valid_aJaI = ...
    lin2accum_aJ > 0 & lin2accum_aJ < (size(grdmag,2) + 1) & ...
    lin2accum_aI > 0 & lin2accum_aI < (size(grdmag,1) + 1);

mask_valid_aJaI_reverse = ~ mask_valid_aJaI;
lin2accum_aJ = lin2accum_aJ .* mask_valid_aJaI + mask_valid_aJaI_reverse;
lin2accum_aI = lin2accum_aI .* mask_valid_aJaI + mask_valid_aJaI_reverse;
clear mask_valid_aJaI_reverse;

% Linear indices (of the votings) into the accumulation array
lin2accum = sub2ind( size(grdmag), lin2accum_aI, lin2accum_aJ );

lin2accum_size = size( lin2accum );
lin2accum = reshape( lin2accum, [numel(lin2accum),1] );
clear lin2accum_aI lin2accum_aJ;

% Weights of the votings, currently using the gradient maginitudes
% but in fact any scheme can be used (application dependent)
weight4accum = ...
    repmat( double(grdmag(grdmasklin)) , [lin2accum_size(2),1] ) .* ...
    mask_valid_aJaI(:);
clear mask_valid_aJaI;

% Build the accumulation array using Matlab function 'accumarray'
accum = accumarray( lin2accum , weight4accum );
accum = [ accum ; zeros( numel(grdmag) - numel(accum) , 1 ) ];
accum = reshape( accum, size(grdmag) );


%%%%%%%% Locating local maxima in the accumulation array %%%%%%%%%%%%

% Stop if no need to locate the center positions of circles
if ~func_compu_cen,
    return;
end
clear lin2accum weight4accum;

% Parameters to locate the local maxima in the accumulation array
% -- Segmentation of 'accum' before locating LM
prm_useaoi = true;
prm_aoithres_s = 2;
prm_aoiminsize = floor(min([ min(size(accum)) * 0.25, ...
    prm_r_range(2) * 1.5 ]));

% -- Filter for searching for local maxima
prm_fltrLM_s = 1.35;
prm_fltrLM_r = ceil( prm_fltrLM_R * 0.6 );
prm_fltrLM_npix = max([ 6, ceil((prm_fltrLM_R/2)^1.8) ]);

% -- Lower bound of the intensity of local maxima
prm_LM_LoBndRa = 0.2;  % minimum ratio of LM to the max of 'accum'

% Smooth the accumulation array
fltr4accum = fltr4accum / sum(fltr4accum(:));
accum = filter2( fltr4accum, accum );

% Select a number of Areas-Of-Interest from the accumulation array
if prm_useaoi,
    % Threshold value for 'accum'
    prm_llm_thres1 = prm_grdthres * prm_aoithres_s;

    % Thresholding over the accumulation array
    accummask = ( accum > prm_llm_thres1 );

    % Segmentation over the mask
    [accumlabel, accum_nRgn] = bwlabel( accummask, 8 );

    % Select AOIs from segmented regions
    accumAOI = ones(0,4);
    for k = 1 : accum_nRgn,
        accumrgn_lin = find( accumlabel == k );
        [accumrgn_IdxI, accumrgn_IdxJ] = ...
            ind2sub( size(accumlabel), accumrgn_lin );
        rgn_top = min( accumrgn_IdxI );
        rgn_bottom = max( accumrgn_IdxI );
        rgn_left = min( accumrgn_IdxJ );
        rgn_right = max( accumrgn_IdxJ );
        % The AOIs selected must satisfy a minimum size
        if ( (rgn_right - rgn_left + 1) >= prm_aoiminsize && ...
                (rgn_bottom - rgn_top + 1) >= prm_aoiminsize ),
            accumAOI = [ accumAOI; ...
                rgn_top, rgn_bottom, rgn_left, rgn_right ];
        end
    end
else
    % Whole accumulation array as the one AOI
    accumAOI = [1, size(accum,1), 1, size(accum,2)];
end

% Thresholding of 'accum' by a lower bound
prm_LM_LoBnd = max(accum(:)) * prm_LM_LoBndRa;

% Build the filter for searching for local maxima
fltr4LM = zeros(2 * prm_fltrLM_R + 1);

[mesh4fLM_x, mesh4fLM_y] = meshgrid(-prm_fltrLM_R : prm_fltrLM_R);
mesh4fLM_r = sqrt( mesh4fLM_x.^2 + mesh4fLM_y.^2 );
fltr4LM_mask = ...
    ( mesh4fLM_r > prm_fltrLM_r & mesh4fLM_r <= prm_fltrLM_R );
fltr4LM = fltr4LM - ...
    fltr4LM_mask * (prm_fltrLM_s / sum(fltr4LM_mask(:)));

if prm_fltrLM_R >= 4,
    fltr4LM_mask = ( mesh4fLM_r < (prm_fltrLM_r - 1) );
else
    fltr4LM_mask = ( mesh4fLM_r < prm_fltrLM_r );
end
fltr4LM = fltr4LM + fltr4LM_mask / sum(fltr4LM_mask(:));

% **** Debug code (begin)
if dbg_on,
    dbg_LMmask = zeros(size(accum));
end
% **** Debug code (end)

% For each of the AOIs selected, locate the local maxima
circen = zeros(0,2);
for k = 1 : size(accumAOI, 1),
    aoi = accumAOI(k,:);    % just for referencing convenience

    % Thresholding of 'accum' by a lower bound
    accumaoi_LBMask = ...
        ( accum(aoi(1):aoi(2), aoi(3):aoi(4)) > prm_LM_LoBnd );

    % Apply the local maxima filter
    candLM = conv2( accum(aoi(1):aoi(2), aoi(3):aoi(4)) , ...
        fltr4LM , 'same' );
    candLM_mask = ( candLM > 0 );

    
    % Clear the margins of 'candLM_mask'
    %Condition added by SPelet
    if size(candLM_mask,1)>prm_fltrLM_R+1 && size(candLM_mask,2)>prm_fltrLM_R+1
        candLM_mask([1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end], :) = 0;
        candLM_mask(:, [1:prm_fltrLM_R, (end-prm_fltrLM_R+1):end]) = 0;
    else
        candLM_mask = 0;
    end

    % **** Debug code (begin)
    if dbg_on,
        dbg_LMmask(aoi(1):aoi(2), aoi(3):aoi(4)) = ...
            dbg_LMmask(aoi(1):aoi(2), aoi(3):aoi(4)) + ...
            accumaoi_LBMask + 2 * candLM_mask;
    end
    % **** Debug code (end)

    % Group the local maxima candidates by adjacency, compute the
    % centroid position for each group and take that as the center
    % of one circle detected
    [candLM_label, candLM_nRgn] = bwlabel( candLM_mask, 8 );

    for ilabel = 1 : candLM_nRgn,
        % Indices (to current AOI) of the pixels in the group
        candgrp_masklin = find( candLM_label == ilabel );
        [candgrp_IdxI, candgrp_IdxJ] = ...
            ind2sub( size(candLM_label) , candgrp_masklin );

        % Indices (to 'accum') of the pixels in the group
        candgrp_IdxI = candgrp_IdxI + ( aoi(1) - 1 );
        candgrp_IdxJ = candgrp_IdxJ + ( aoi(3) - 1 );
        candgrp_idx2acm = ...
            sub2ind( size(accum) , candgrp_IdxI , candgrp_IdxJ );

        % Minimum number of qulified pixels in the group
        if sum(accumaoi_LBMask(candgrp_masklin)) < prm_fltrLM_npix,
            continue;
        end

        % Compute the centroid position
        candgrp_acmsum = sum( accum(candgrp_idx2acm) );
        cc_x = sum( candgrp_IdxJ .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
        cc_y = sum( candgrp_IdxI .* accum(candgrp_idx2acm) ) / ...
            candgrp_acmsum;
        circen = [circen; cc_x, cc_y];
    end
end

% **** Debug code (begin)
if dbg_on,
    figure(dbg_bfigno); imagesc(dbg_LMmask); axis image;
    title('Generated map of local maxima');
    if size(accumAOI, 1) == 1,
        figure(dbg_bfigno+1);
        surf(candLM, 'EdgeColor', 'none'); axis ij;
        title('Accumulation array after local maximum filtering');
    end
end
% **** Debug code (end)


%%%%%%%% Estimation of the Radii of Circles %%%%%%%%%%%%

% Stop if no need to estimate the radii of circles
if ~func_compu_radii,
    varargout{1} = circen;
    return;
end

% Parameters for the estimation of the radii of circles
fltr4SgnCv = [2 1 1];
fltr4SgnCv = fltr4SgnCv / sum(fltr4SgnCv);

% Find circle's radius using its signature curve
cirrad = zeros( size(circen,1), 1 );

for k = 1 : size(circen,1),
    % Neighborhood region of the circle for building the sgn. curve
    circen_round = round( circen(k,:) );
    SCvR_I0 = circen_round(2) - prm_r_range(2) - 1;
    if SCvR_I0 < 1,
        SCvR_I0 = 1;
    end
    SCvR_I1 = circen_round(2) + prm_r_range(2) + 1;
    if SCvR_I1 > size(grdx,1),
        SCvR_I1 = size(grdx,1);
    end
    SCvR_J0 = circen_round(1) - prm_r_range(2) - 1;
    if SCvR_J0 < 1,
        SCvR_J0 = 1;
    end
    SCvR_J1 = circen_round(1) + prm_r_range(2) + 1;
    if SCvR_J1 > size(grdx,2),
        SCvR_J1 = size(grdx,2);
    end

    % Build the sgn. curve
    SgnCvMat_dx = repmat( (SCvR_J0:SCvR_J1) - circen(k,1) , ...
        [SCvR_I1 - SCvR_I0 + 1 , 1] );
    SgnCvMat_dy = repmat( (SCvR_I0:SCvR_I1)' - circen(k,2) , ...
        [1 , SCvR_J1 - SCvR_J0 + 1] );
    SgnCvMat_r = sqrt( SgnCvMat_dx .^2 + SgnCvMat_dy .^2 );
    SgnCvMat_rp1 = round(SgnCvMat_r) + 1;

    f4SgnCv = abs( ...
        double(grdx(SCvR_I0:SCvR_I1, SCvR_J0:SCvR_J1)) .* SgnCvMat_dx + ...
        double(grdy(SCvR_I0:SCvR_I1, SCvR_J0:SCvR_J1)) .* SgnCvMat_dy ...
        ) ./ SgnCvMat_r;
    SgnCv = accumarray( SgnCvMat_rp1(:) , f4SgnCv(:) );

    SgnCv_Cnt = accumarray( SgnCvMat_rp1(:) , ones(numel(f4SgnCv),1) );
    SgnCv_Cnt = SgnCv_Cnt + (SgnCv_Cnt == 0);
    SgnCv = SgnCv ./ SgnCv_Cnt;

    % Suppress the undesired entries in the sgn. curve
    % -- Radii that correspond to short arcs
    SgnCv = SgnCv .* ( SgnCv_Cnt >= (pi/4 * [0:(numel(SgnCv_Cnt)-1)]') );
    % -- Radii that are out of the given range
    SgnCv( 1 : (round(prm_r_range(1))+1) ) = 0;
    SgnCv( (round(prm_r_range(2))+1) : end ) = 0;

    % Get rid of the zero radius entry in the array
    SgnCv = SgnCv(2:end);
    % Smooth the sgn. curve
    SgnCv = filtfilt( fltr4SgnCv , [1] , SgnCv );

    % Get the maximum value in the sgn. curve
    SgnCv_max = max(SgnCv);
    if SgnCv_max <= 0,
        cirrad(k) = 0;
        continue;
    end

    % Find the local maxima in sgn. curve by 1st order derivatives
    % -- Mark the ascending edges in the sgn. curve as 1s and
    % -- descending edges as 0s
    SgnCv_AscEdg = ( SgnCv(2:end) - SgnCv(1:(end-1)) ) > 0;
    % -- Mark the transition (ascending to descending) regions
    SgnCv_LMmask = [ 0; 0; SgnCv_AscEdg(1:(end-2)) ] & (~SgnCv_AscEdg);
    SgnCv_LMmask = SgnCv_LMmask & [ SgnCv_LMmask(2:end) ; 0 ];

    % Incorporate the minimum value requirement
    SgnCv_LMmask = SgnCv_LMmask & ...
        ( SgnCv(1:(end-1)) >= (prm_multirad * SgnCv_max) );
    % Get the positions of the peaks
    SgnCv_LMPos = sort( find(SgnCv_LMmask) );

    % Save the detected radii
    if isempty(SgnCv_LMPos),
        cirrad(k) = 0;
    else
        cirrad(k) = SgnCv_LMPos(end);
        for i_radii = (length(SgnCv_LMPos) - 1) : -1 : 1,
            circen = [ circen; circen(k,:) ];
            cirrad = [ cirrad; SgnCv_LMPos(i_radii) ];
        end
    end
end

% Output
varargout{1} = circen;
varargout{2} = cirrad;
if nargout > 3,
    varargout{3} = dbg_LMmask;
end

%%

function DrawCircle(x, y, r, nseg, S)
% Draw a circle on the current figure using ploylines
%
%  DrawCircle(x, y, r, nseg, S)
%  A simple function for drawing a circle on graph.
%
%  INPUT: (x, y, r, nseg, S)
%  x, y:    Center of the circle
%  r:       Radius of the circle
%  nseg:    Number of segments for the circle
%  S:       Colors, plot symbols and line types
%
%  OUTPUT: None
%
%  BUG REPORT:
%  Please send your bug reports, comments and suggestions to
%  pengtao@glue.umd.edu . Thanks.

%  Author:  Tao Peng
%           Department of Mechanical Engineering
%           University of Maryland, College Park, Maryland 20742, USA
%           pengtao@glue.umd.edu
%  Version: alpha       Revision: Jan. 10, 2006


theta = 0 : (2 * pi / nseg) : (2 * pi);
pline_x = r * cos(theta) + x;
pline_y = r * sin(theta) + y;

plot(pline_x, pline_y, S);

%%

function ImOut = CutWithBorder(ImIn, Corners, Border);
%corner: [Ymin, Ymax, Xmin, Xmax]
ImMid = ImIn(Corners(1):Corners(2),Corners(3):Corners(4));

ImOut = zeros(size(ImMid,1)+2*Border,size(ImMid,2)+2*Border);

ImOut(Border+1:size(ImMid,1)+Border,Border+1:size(ImMid,2)+Border)  = ImMid;


