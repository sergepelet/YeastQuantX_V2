function Var = SegmentNuclBF(Var, CallNum)
tic
Debug = 0;     %Set to 1 to display segmentation images and 0 not to

if Debug == 1
    warning on
else
    warning off
end

%load file containing segmentation parameters
SegParaFullName = [Var.Analysis.SegAroundPara{CallNum}, '_', num2str(Var.Experiment.Objective),'x_bin', num2str(Var.Experiment.Bin)];

ParaNum = [];
%Identify correct segmentation paramters
for i = 1:length(Var.SegmentationParameters)
    if strcmp(SegParaFullName, Var.SegmentationParameters(i).FullName)
        ParaNum = i;
    end
end



%min and Max diameter of objects
MinDiam = Var.SegmentationParameters(ParaNum).MinDiameter;
MaxDiam = Var.SegmentationParameters(ParaNum).MaxDiameter;
IncludeEdge = Var.SegmentationParameters(ParaNum).Touch;
NuclCytoPix = Var.SegmentationParameters(ParaNum).NuclCytoDistance;

%Parameter to increase the size of the Bounding box
GrowBoundingBox = Var.SegmentationParameters(ParaNum).BoundingBox; %30;

%Parameter for object dilation
SmallDilate = Var.SegmentationParameters(ParaNum).SmallDilation; %3;
LargeDilate = Var.SegmentationParameters(ParaNum).LargeDilation; %6
%Load Nucl Object and BF images
NuclImg = Var.Img.(Var.Analysis.SegYFImg{CallNum});
NuclObj = Var.Img.(Var.Analysis.SegYFImgOut{CallNum});

InFocus = Var.Img.(Var.Analysis.SegPhiImg{CallNum});
InFocus = (InFocus-min(InFocus(:)))./(max(InFocus(:))-min(InFocus(:)));
% Load and enhance BF out of focus Image
ZOffImg = Var.Img.(Var.Analysis.SegPhiSecImg{CallNum});
ZOffImg = (ZOffImg-min(ZOffImg(:)))./(max(ZOffImg(:))-min(ZOffImg(:)));

%Get Nuclei Labels
NuclLabel = Var.Measurements.(Var.Analysis.SegYFImgOut{CallNum}).SegLabel;



%% Use NuclImg and BF image to find cell boundaries
%Label Nucl image
NuclObjStuct = bwconncomp(NuclObj);
NuclObj = labelmatrix(NuclObjStuct);
NumNucl = NuclObjStuct.NumObjects;
StartNucl = NumNucl;
NuclProps = regionprops(NuclObj,'BoundingBox','ConvexArea','ConvexImage', 'Area');


%Border parameter sets the number of pixels added on each sides of the
%bounding box defined around the large areas
Border = 0;


CombCell  = zeros(size(NuclObj));
CombNucl = zeros(size(NuclObj));
CombOFFCell = zeros(size(NuclObj));

if NumNucl > 0
    for N = 1:NumNucl
        %N = N
        % Create image containing only one large area
        %use the bounding box to define the corners of the small image
        ULCorner = round(NuclProps(N).BoundingBox(1:2));
        Width = NuclProps(N).BoundingBox(3:4);
        Corners = [ULCorner(2)-GrowBoundingBox,ULCorner(2)+Width(2)+GrowBoundingBox,...
            ULCorner(1)-GrowBoundingBox,ULCorner(1)+Width(1)+GrowBoundingBox ];
        %Check if corners are within the size of the full image
        if Corners(3) < 1
            Corners(3) = 1;
        end
        if Corners(4) > size(NuclObj,2)
            Corners(4) = size(NuclObj,2);
        end
        if Corners(1) < 1
            Corners(1) = 1;
        end
        if Corners(2) > size(NuclObj,1)
            Corners(2) = size(NuclObj,1);
        end
        
        %Cut full images to keep only the small region that contain the
        %area of interest
        %And Remove Pixels belonging to other objects
        CutNucl = CutWithBorder(NuclObj, Corners, Border);
        CutNucl(CutNucl~= N) = 0;
        CutNucl(CutNucl>0) = 1;
        %assignin('base',['CutNucl', num2str(N)], CutNucl)
        %Cut BF image
        CutZOff = CutWithBorder(ZOffImg, Corners, Border);
        CutFocus = CutWithBorder(InFocus, Corners, Border);
        % assignin('base',['CutZOff_', num2str(N)], CutZOff)
        %Normalize image and apply levelset
        IntLimit = 0.02;
        CutNorm = imadjust(CutZOff,stretchlim(CutZOff, [IntLimit,1-IntLimit]),[]);
        % CutFilter = wiener2(CutNorm,[10 10]);
        %CutNorm = CutZOff/max(CutZOff(:));
        CutSeg = region_seg(CutNorm, CutNucl, 250,0.2);
        % CutSegFilter = region_seg(CutFilter, CutNucl, 250,0.2);
        if Debug
            figure(3000+N)
            % CutSegFilter = double(CutSegFilter);
            % CutSegFilter(CutSeg)= 3;
            subplot(2,2,1); imagesc(CutNucl), title(['CutNucl ', num2str(N)])
            
            subplot(2,2,2); imagesc(CutZOff), title('CutZOff')
            subplot(2,2,3); imagesc(CutFocus), title('CutFocus')
            subplot(2,2,4); imagesc(CutSeg), title('CutSeg')
            pause(1)
        end
        
        % PROBLEM IF NO AREA FOUND FOR NOW SET TO NUCL AREA
        if max(CutSeg(:)) ==0
            
            
            %Try to change the normalization limit for BF image
            IntLimit = IntLimit/2;
            %Redo level set
            CutNorm = imadjust(CutZOff,stretchlim(CutZOff, [IntLimit,1-IntLimit]),[]);
            CutSeg = region_seg(CutNorm, CutNucl, 250,0.2);
            if Debug
                figure(300+N)
                subplot(2,2,1); imagesc(CutNucl), title('CutNucl')
                
                subplot(2,2,2); imagesc(CutZOff), title('CutZOff')
                subplot(2,2,3); imagesc(CutSeg), title('CutSeg')
                subplot(2,2,4); imagesc(CutFocus), title('CutFocus')
            end
            if max(CutSeg(:)) ==0
                CutSeg = CutNucl;
                warning(['No area found for Nucl Num: ', num2str(N)])
                if Debug
                    figure(300+N)
                    subplot(2,2,1); imagesc(CutNucl), title('CutNucl')
                    
                    subplot(2,2,2); imagesc(CutZOff), title('CutZOff')
                    subplot(2,2,3); imagesc(CutSeg), title('CutSeg')
                    subplot(2,2,4); imagesc(CutFocus), title('CutFocus')
                    pause(1)
                    
                    assignin('base',['CutNucl_', num2str(N)], CutNucl)
                    assignin('base',['CutZOff_', num2str(N)], CutZOff)
                end
            end
            
        end
        
        %         assignin('base',['CutNucl_', num2str(N)], CutNucl)
        %                 assignin('base',['CutZOff_', num2str(N)], CutZOff)
        %                 error('BLA')
        %Set Pixels belonging to the nucleus to 1 in CutSeg Image
        CutSeg(CutNucl>0) = 1;
        
        
        %Get convex image from segmented object
        CutProps = regionprops(CutSeg, 'ConvexImage', 'BoundingBox');
        %Generate CutConv image based on the bounding box of the object
        CutConv = zeros(size(CutNucl));
        CutConv(round(CutProps(1).BoundingBox(2)):round(CutProps(1).BoundingBox(2))+CutProps(1).BoundingBox(4)-1,...
            round(CutProps(1).BoundingBox(1)):round(CutProps(1).BoundingBox(1))+CutProps(1).BoundingBox(3)-1) = CutProps(1).ConvexImage;
        
        
        
        
        %Dilate Image
        SE = strel('disk',SmallDilate);
        CutConvDil = imdilate(CutConv,SE);
        %Normalize Image
        IntLimit = 0.05;
        CutFocusNorm = imadjust(CutFocus,stretchlim(CutFocus, [IntLimit,1-IntLimit]),[]);
        %Use level-set to expand CutConv Image based on Image in focus
        CutSegFocus = region_seg(CutFocusNorm, CutConvDil, 250,0.2);
        
        %Dilate Image with larger radius
        SE = strel('disk',LargeDilate);
        CutConvDilHi = imdilate(CutConv,SE);
        
        %Exclude pixels which lie outside
        CutSegFocus(CutConvDilHi==0) = 0;
        CutSegFocus(CutConv>0) = 1;
        
        %Get convex image and transfer
        CutProps = regionprops(CutSegFocus, 'ConvexImage', 'BoundingBox');
        CutConvFocus = zeros(size(CutNucl));
        CutConvFocus(round(CutProps(1).BoundingBox(2)):round(CutProps(1).BoundingBox(2))+CutProps(1).BoundingBox(4)-1,...
            round(CutProps(1).BoundingBox(1)):round(CutProps(1).BoundingBox(1))+CutProps(1).BoundingBox(3)-1) =CutProps(1).ConvexImage;
        
        %Smooth edges by opening
        SE = strel('disk',2);
        CutConvFocus = imopen(CutConvFocus,SE);
        
        %         figure(7000 +N)
        %                  subplot(2,2,1); imagesc(CutConv)
        %
        %                  subplot(2,2,2); imagesc(CutConvFocus)
        %                  pause(1)
        %test overlap between Nucleus and cell
        %Get Nucleus pixels
        NuclPix = find(CutNucl > 0);
        if min(CutConvFocus(NuclPix)) == 1;     %Test if all pixels are one
            %ALL GOOD Whole Nucl is in the cell
            %Get Cell pixels
            CellPix = find(CutConvFocus > 0);
            % fprintf(['NumPix = ' num2str(length(NuclPix)), 'NumCell = ' num2str(length(CellPix)), 'Ratio = ', num2str(length(NuclPix)/length(CellPix)),'\n'])
            %Test if cell is bigger than nucleus
            if length(NuclPix)/length(CellPix) > 0.5
                %PB cell is too small
                %Remove cell from analysis
                warning(['Cell to small for Nulc: ', num2str(N), ' NumPix = ' num2str(length(NuclPix)), ' NumCell = ' num2str(length(CellPix)), ' Ratio = ', num2str(length(NuclPix)/length(CellPix))])
                continue
            end
            
        elseif min(CutConvFocus(NuclPix)) == 0;     %Test if some nuclear Pix are zero (outside of object)
            %PB: some or all pix of Nucl are outside of the cell
            %Get overlapping pixels
            OverlapPix = find(CutConvFocus(NuclPix) >0);
            NoOverlap = find(CutConvFocus(NuclPix) == 0);
            %Test if Overlap is more than half of the nucl
            if length(OverlapPix) >= length(NoOverlap)
                % if Yes => just remove the non-overlapping pixels
                %                  figure(5000 +N)
                %                  subplot(2,2,1); imagesc(CutNucl)
                %                  subplot(2,2,3); imagesc(CutNucl+CutConvFocus)
                CutNucl(NoOverlap) = 0;
                %       subplot(2,2,2); imagesc(CutNucl)
                
                %       pause(1)
            else
                % If no or small overlap between Nucl and Focused object
                %Try to use the out of focus oobject
                %Test if all Nucl pixel belong to the object
                if min(CutConv(NuclPix)) == 1
                    % Test if the cut conv object is at least twice as
                    % large as the nucl
                    if length(find(CutNucl == 1))*2 <= length(find(CutConv == 1))
                        %Set object to out of focus object which is the
                        %best guess for the cell shape
                        CutConvFocus = CutConv;
                    else
                        %If not the case remove this object
                        warning(['Bad Overlap between Nucl and cell for Nulc: ', num2str(N),' CutConv too small'])
                        continue
                    end
                else
                %    If not the case remove this object
                    if Debug
                        figure(800+N)
                        
                        Foc_Nucl = CutConvFocus;
                        Foc_Nucl(CutNucl ==1) = 2;
                        Off_Nucl = CutConv;
                        Off_Nucl(CutNucl ==1) = 2;
                        subplot(2,2,1); imagesc(Off_Nucl), title('Off_Nucl')
                        subplot(2,2,3); imagesc(Foc_Nucl), title('Foc_Nucl')
                        P_Conv = bwperim(CutConv);
                        P_CutFocusNorm = CutFocusNorm;
                        P_CutFocusNorm(P_Conv) = 1;
                        
                        subplot(2,2,2); imagesc(P_CutFocusNorm), title('P_CutFocusNorm_OFF')
                        P_Foc = bwperim(CutConvFocus);
                        P_CutFocusNorm = CutFocusNorm;
                        P_CutFocusNorm(P_Foc) = 1;
                        subplot(2,2,4); imagesc(P_CutFocusNorm), title('P_CutFocusNorm focus')
                        pause(1)
                    end
                    warning(['Bad Overlap between Nucl and cell for Nulc: ', num2str(N),' Nucl outside of CutConv'])
                    continue
                end
                
                
                
                
                
                if Debug
                    figure(400+N)
                    subplot(2,2,1); imagesc(CutNucl), title('CutNucl')
                    
                    subplot(2,2,2); imagesc(CutFocusNorm), title('CutFocusNorm')
                    subplot(2,2,3); imagesc(CutConv), title('CutConv')
                    subplot(2,2,4); imagesc(CutConvFocus), title('CutConvFocus')
                    
                    assignin('base',['CutNucl_', num2str(N)], CutNucl)
                    assignin('base',['CutConv_', num2str(N)], CutConv)
                    assignin('base',['CutFocus_', num2str(N)], CutFocus)
                    assignin('base',['CutConvDil_', num2str(N)], CutConvDil)
                    
                    figure(800+N)
                    
                    Foc_Nucl = CutConvFocus;
                    Foc_Nucl(CutNucl ==1) = 2;
                    Off_Nucl = CutConv;
                    Off_Nucl(CutNucl ==1) = 2;
                    subplot(2,2,1); imagesc(Off_Nucl), title('Off_Nucl')
                    subplot(2,2,3); imagesc(Foc_Nucl), title('Foc_Nucl')
                    P_Conv = bwperim(CutConv);
                    P_CutFocusNorm = CutFocusNorm;
                    P_CutFocusNorm(P_Conv) = 1;
                    
                    subplot(2,2,2); imagesc(P_CutFocusNorm), title('P_CutFocusNorm_OFF')
                    P_Foc = bwperim(CutConvFocus);
                    P_CutFocusNorm = CutFocusNorm;
                    P_CutFocusNorm(P_Foc) = 1;
                    subplot(2,2,4); imagesc(P_CutFocusNorm), title('P_CutFocusNorm focus')
                    pause(1)
                    
                end
                warning(['Bad Overlap between Nucl and cell for Nulc: ', num2str(N)])
                continue
            end
        end
        
        %Test size of Cell compared to Nucl
        
        %See if it is large enough to generate the  cyto object
        %Dilate image by NuclCytoPix distance
        SE = strel('disk',NuclCytoPix);
        CutNuclGrow = imdilate(CutNucl,SE);
        CutCyto = CutConvFocus;
        CutCyto(CutNuclGrow==1) = 0;
        %CytopNulc = length(find(CutCyto == 1)) + length(find(CutNucl == 1))/1000
        if length(find(CutNucl == 1)) > length(find(CutCyto == 1))
            CutCyto(CutNuclGrow==1) = 2;
            if Debug
                figure(700+N)
                subplot(2,1,1); imagesc(CutCyto)
                subplot(2,1,2); imagesc(CutFocusNorm), title('CutFocusNorm')
            end
            warning(['Cyto area smaller than Nucl area for Nucl: ', num2str(N)])
            continue
        end
        
        
        %TRANSFER CUT IMAGE to the full-size image
        %Get indices of the Nucleus
        SegInd = find(CutNucl > 0);
        %Convert indices to subscripts
        [SegI,SegJ] = ind2sub(size(CutNucl),SegInd);
        
        %Find subscript that would fall outside of the border of the image
        LowI = find(SegI+Corners(1)-Border-1<1);
        HiI = find(SegI+Corners(1)-Border-1>size(CombNucl,1));
        LowJ = find(SegJ+Corners(3)-Border-1<1);
        HiJ = find(SegJ+Corners(3)-Border-1>size(CombNucl,2));
        
        %Collect the indices of all outlier subscripts found
        AllOUT = [];
        if ~isempty(LowI)
            AllOUT = [AllOUT;LowI];
        end
        if ~isempty(HiI)
            AllOUT = [AllOUT;HiI];
        end
        if ~isempty(LowJ)
            AllOUT = [AllOUT;LowJ];
        end
        if ~isempty(HiJ)
            AllOUT = [AllOUT;HiJ];
        end
        %Find indices that are not present in outliers
        AllIN = setdiff([1:length(SegI)], AllOUT);
        SegI = SegI(AllIN);
        SegJ = SegJ(AllIN);
        %Convert back to indices for the full-size image.
        SegFullInd = sub2ind(size(CombNucl),SegI+Corners(1)-Border-1,SegJ+Corners(3)-Border-1);   %corner: [Ymin, Ymax, Xmin, Xmax]
        %Set the value of the pixels in final image to N
        CombNucl(SegFullInd) = N;
        
        
        %Get indices of the Cell
        SegInd = find(CutConvFocus > 0);
        %Convert indices to subscripts
        [SegI,SegJ] = ind2sub(size(CutConvFocus),SegInd);
        
        %Find subscript that would fall outside of the border of the image
        LowI = find(SegI+Corners(1)-Border-1<1);
        HiI = find(SegI+Corners(1)-Border-1>size(CombCell,1));
        LowJ = find(SegJ+Corners(3)-Border-1<1);
        HiJ = find(SegJ+Corners(3)-Border-1>size(CombCell,2));
        %Collect the indices of all outlier subscripts found
        AllOUT = [];
        if ~isempty(LowI)
            AllOUT = [AllOUT;LowI];
        end
        if ~isempty(HiI)
            AllOUT = [AllOUT;HiI];
        end
        if ~isempty(LowJ)
            AllOUT = [AllOUT;LowJ];
        end
        if ~isempty(HiJ)
            AllOUT = [AllOUT;HiJ];
        end
        
        %Find indices that are not present in outliers
        AllIN = setdiff([1:length(SegI)], AllOUT);
        SegI = SegI(AllIN);
        SegJ = SegJ(AllIN);
        %Convert back to indices for the full-size image.
        SegFullInd = sub2ind(size(CombCell),SegI+Corners(1)-Border-1,SegJ+Corners(3)-Border-1);   %corner: [Ymin, Ymax, Xmin, Xmax]
        %Set the value of the pixels in final image to N
        CombCell(SegFullInd) = N;
        
        
        %Get indices of the Cell
        SegInd = find(CutConv > 0);
        %Convert indices to subscripts
        [SegI,SegJ] = ind2sub(size(CutConv),SegInd);
        
        %Find subscript that would fall outside of the Border of the image
        LowI = find(SegI+Corners(1)-Border-1<1);
        HiI = find(SegI+Corners(1)-Border-1>size(CombOFFCell,1));
        LowJ = find(SegJ+Corners(3)-Border-1<1);
        HiJ = find(SegJ+Corners(3)-Border-1>size(CombOFFCell,2));
        %Collect the indices of all outlier subscripts found
        AllOUT = [];
        if ~isempty(LowI)
            AllOUT = [AllOUT;LowI];
        end
        if ~isempty(HiI)
            AllOUT = [AllOUT;HiI];
        end
        if ~isempty(LowJ)
            AllOUT = [AllOUT;LowJ];
        end
        if ~isempty(HiJ)
            AllOUT = [AllOUT;HiJ];
        end
        
        %Find indices that are not present in outliers
        AllIN = setdiff([1:length(SegI)], AllOUT);
        SegI = SegI(AllIN);
        SegJ = SegJ(AllIN);
        %Convert back to indices for the full-size image.
        SegFullInd = sub2ind(size(CombOFFCell),SegI+Corners(1)-Border-1,SegJ+Corners(3)-Border-1);   %corner: [Ymin, Ymax, Xmin, Xmax]
        %Set the value of the pixels in final image to N
        CombOFFCell(SegFullInd) = N;
        
    end
    %     if Debug
    %         figure(200)
    %         imagesc(CombCell)
    %                 figure(203)
    %         imagesc(CombOFFCell)
    %     end
    
    %Label nuclei
    NuclStruct = bwconncomp(CombNucl);
    LabelNucl = double(labelmatrix(NuclStruct));
    NewNumNucl= NuclStruct.NumObjects;
    NuclProps = regionprops(LabelNucl, 'Centroid');
    %Label cells
    CellStruct = bwconncomp(CombCell);
    LabelCell = double(labelmatrix(CellStruct));
    NumCell = CellStruct.NumObjects;
    CellProps = regionprops(LabelCell, 'BoundingBox', 'Image', 'PixelIdxList');
    
    %Check overlap between different cells
    %OverlapCell = zeros(size(CombCell));
    NuclConflict = [];
    CellConflict = [];
    for C = 1:NumCell
        %Get Number of nucleus in object
        NuclID = unique(LabelNucl(CellStruct.PixelIdxList{C}));
        NuclID = NuclID(NuclID>0);
        NuclFound = length(NuclID);
        switch NuclFound
            case 0  % No overlap with nucleus Remove cell
                warning(['Error! no overlap with nucleus for Cell:', num2str(C)])
                %Set  Cell to Zero
                
            case 1 %ALL good
                %CombCell(CellStruct.PixelIdxList{C}) = NuclID;
            otherwise   % More than one nucl associated with region need to split
                %Add cells to overlap Cell
                %OverlapCell(CellStruct.PixelIdxList{C}) = 1;
                NuclConflict = [NuclConflict, NuclID'];
                CellConflict = [CellConflict, C.*ones(1,NuclFound)];
        end
    end
    
    % CombCellDup = CombCell;
    %Go back to every cell generating a conflict
    if ~isempty(CellConflict)
        %Get list of all conflict cells
        CellList = unique(CellConflict);
        %Loop through all Conflict cells
        for C = CellList
            %C = C
            IndC = find(CellConflict == C);
            NumNucl = length(IndC);
            NuclID = NuclConflict(IndC);
            
            
            GroupCellImg = CellProps(C).Image;
            
            %Cut Nucl Image to extract nuclei
            C_ULCorner = round(CellProps(C).BoundingBox(1:2));
            C_Width = CellProps(C).BoundingBox(3:4);
            C_Corners = [C_ULCorner(2),C_ULCorner(2)+C_Width(2)-1,...
                C_ULCorner(1),C_ULCorner(1)+C_Width(1)-1 ];
            
            C_CutNuclImg = CutWithBorder(LabelNucl, C_Corners, 0);
            
            %C_CutCellImg = CutWithBorder(CombCell, C_Corners, 0);
            
            %Create small Img with nucl set to one
            ConflictNuclImg = zeros(size(C_CutNuclImg));
            for NID = NuclID
                ConflictNuclImg(C_CutNuclImg==NID) = 1;
            end
            
            %                                 figure(600+C)
            %                     subplot(2,2,1); imagesc(ConflictNuclImg), title(['ConflictNuclImg', num2str(NuclID(1)), '  ', num2str(NuclID(2))])
            %
            %                     subplot(2,2,2); imagesc(GroupCellImg), title('GroupCellImg')
            
            
            
            %%Use watershed simple distance watershed to split cells containing multiple nuclei
            %Calculate distance transform
            DistCell = bwdist(~GroupCellImg);
            DistCell = -DistCell;
            %DistCell(~GroupCellImg) = -Inf;
            %Some smoothing to avoid too many local minima
            DistCellMin = imhmin(DistCell,1);
            WaterCell = watershed(DistCellMin);
            
            %Apply the splitting from Watershed
            SplitCell = double(GroupCellImg);
            SplitCell(WaterCell==0) = 0;
            %  subplot(2,2,3); imagesc(SplitCell), title('SplitCell')
            %Remove area not overlapping with nucleus
            SplitCell(ConflictNuclImg>0) = 2;
            % subplot(2,2,3); imagesc(SplitCell), title('SplitCell')
            SplitCell = imextendedmax(SplitCell,1);
            % subplot(2,2,4); imagesc(SplitCell), title('SplitCell Imextend')
            %Define regions in Split cellImg
            SplitStruct = bwconncomp(SplitCell);
            NumSplitCells = SplitStruct.NumObjects;
            
            %  subplot(2,2,3); imagesc(SplitCell), title('SplitCell after imextendmax')
            % Check if Number of cell and nucleus match
            if NumSplitCells == NumNucl
                SplitInd = find(SplitCell > 0);
                %Convert indices to subscripts
                [SegI,SegJ] = ind2sub(size(SplitCell),SplitInd);
                
                %Convert back to indices for the full-size image.
                SegFullInd = sub2ind(size(LabelCell),SegI+C_Corners(1)-1,SegJ+C_Corners(3)-1);   %corner: [Ymin, Ymax, Xmin, Xmax]
                %Remove Cell from CombCell
                CombCell(CellProps(C).PixelIdxList) =0;
                
                %Set the value of the pixels in cell image to N
                CombCell(SegFullInd) = 1;
                %                CombCellDup(SegFullInd) = CombCellDup(SegFullInd)+5;
                
            else %IF splitting based on cell shape did not work, use Nuclei to split the image
                
                %Calculate distance transform
                DistNucl = bwdist(ConflictNuclImg);
                
                %Some smoothing to avoid too many local minima
                DistNuclMin = imhmin(DistNucl,1);
                WaterNucl = watershed(DistNuclMin);
                if max(WaterNucl(:)) < NumNucl
                    SE = strel('disk',1);
                    ShrinkConflictNuclImg = imerode(ConflictNuclImg, SE);
                    %Calculate distance transform
                    DistNucl = bwdist(ShrinkConflictNuclImg);
                    
                    %Some smoothing to avoid too many local minima
                    DistNuclMin = imhmin(DistNucl,1);
                    WaterNucl = watershed(DistNuclMin);
                end
                
                %                 subplot(2,2,4); imagesc(DistNuclMin), title('DistNuclMin')
                %                 subplot(2,2,3); imagesc(WaterNucl), title('WaterNucl')
                
                
                %Apply the splitting from Watershed
                SplitCell = GroupCellImg;
                SplitCell(WaterNucl==0) = 0;
                
                
                %Remove area not overlapping with nucleus
                SplitCell = double(SplitCell);
                SplitCell(ConflictNuclImg>0) = 2;
                % subplot(2,2,3); imagesc(SplitCell), title('SplitCell + Nucl')
                
                SplitCell = imextendedmax(SplitCell,1);
                
                %Define regions in Split cellImg
                SplitStruct = bwconncomp(SplitCell);
                NumSplitCells = SplitStruct.NumObjects;
                %subplot(2,2,4); imagesc(SplitCell), title('SplitCell Nucl watershed')
                
                if NumSplitCells == NumNucl
                    SplitInd = find(SplitCell > 0);
                    %Convert indices to subscripts
                    [SegI,SegJ] = ind2sub(size(SplitCell),SplitInd);
                    
                    %Convert back to indices for the full-size image.
                    SegFullInd = sub2ind(size(LabelCell),SegI+C_Corners(1)-1,SegJ+C_Corners(3)-1);   %corner: [Ymin, Ymax, Xmin, Xmax]
                    
                    %Remove Cell from CombCell
                    CombCell(CellProps(C).PixelIdxList) =0;
                    
                    %Set the value of the pixels in cell image to N
                    CombCell(SegFullInd) = 1;
                    % CombCellDup(SegFullInd) = 1;
                else
                    %SHOULD get rid of cells that cannot be split correctly
                    %Keep for now for debugging.
                    fprintf('Pb with Nucl Watershed split')
                end
                
            end
        end
    end
    
    
    
    
    
    %Re-label Cell image
    CellStruct = bwconncomp(CombCell);
    LabelCell = labelmatrix(CellStruct);
    NumCell = CellStruct.NumObjects;
    CellProps = regionprops(LabelCell, 'Centroid');
    
    if Debug
        figure(201)
        CellNucl = LabelCell;
        CellNucl(LabelNucl>0) = CellNucl(LabelNucl>0)+5;
        imagesc(CellNucl)
    end
    
    %If mismatch between Nucl and Cell
    if NumCell ~= NewNumNucl
        %Loop through all cells
        NuclIndList = [];
        for O = 1:max(LabelNucl(:))
            %Find Cell Ind belonging to a Nucleus
            CellInd = unique(LabelCell(LabelNucl ==O));
            CellInd = CellInd(CellInd~= 0);
            %If more than one cell belong to given nucl
            if length(CellInd) > 1
                %Remove cell objects
                for I = 1:length(CellInd)
                    LabelCell(LabelCell == CellInd(I)) = 0;
                end
                %Remove Nucl Object
                LabelNucl(LabelNucl ==O) = 0;
            end
            
            %Re-label Cell image
            CellStruct = bwconncomp(LabelCell);
            LabelCell = labelmatrix(CellStruct);
            NumCell = CellStruct.NumObjects;
            CellProps = regionprops(LabelCell, 'Centroid');
            
            %Relabel NuclImg
            NuclStruct = bwconncomp(LabelNucl);
            LabelNucl = labelmatrix(NuclStruct);
            NewNumNucl = NuclStruct.NumObjects;
            NuclProps = regionprops(LabelNucl, 'Centroid');
        end
        
        if NumCell ~= NewNumNucl
            error(['Error number of cells (', num2str(NumCell), ') NOT EQUAL to number of nuclei (', num2str(NewNumNucl), ')'])
        end
    end
    
    %% Remove Object touching the border of the image
    
    if Debug
        figure(600)
        imagesc(LabelCell); title('BEOFRE')
    end
    
    if strcmp(IncludeEdge,'No')
        LabelCell = imclearborder(LabelCell);
    end
    
    
    %Re-label Cell image
    CellStruct = bwconncomp(LabelCell);
    LabelCell = labelmatrix(CellStruct);
    NumCell = CellStruct.NumObjects;
    CellProps = regionprops(LabelCell, 'Centroid');
    if Debug
        figure(601)
        imagesc(LabelCell); title('After')
    end
    
    %Remove Cells touching border from nucleus image
    LabelNucl(LabelCell == 0) = 0;
    NuclStruct = bwconncomp(LabelNucl);
    LabelNucl = labelmatrix(NuclStruct);
    FinalNuclNum = NuclStruct.NumObjects;
    NuclProps = regionprops(LabelNucl, 'Centroid');
    
    
    
    
    %Transfer Nucl Label to Cell image
    FinalCell = zeros(size(CombCell));
    NuclLabelOUT = zeros(FinalNuclNum,1);
    NuclCentr1 = zeros(FinalNuclNum,1);
    NuclCentr2 = zeros(FinalNuclNum,1);
    CellCentr1 = zeros(FinalNuclNum,1);
    CellCentr2 = zeros(FinalNuclNum,1);
    if Debug
        figure(800)
        imagesc(LabelNucl)
        figure(801)
        imagesc(LabelCell)
    end
    for N = 1:FinalNuclNum
        %Get First Nuclear Pixel
        NuclPix = NuclStruct.PixelIdxList{N}(1);
        %Get Nucl and Cell labels at that pixel
        NuclID = LabelNucl(NuclPix);
        CellID = LabelCell(NuclPix);
        if Debug
            [NY,NX] = ind2sub(size(LabelNucl),NuclPix);
            
            %         figure(800)
            %         hold on
            %         plot(NX,NY, 'or')
            %
            %            figure(801)
            %         hold on
            %         plot(NX,NY, 'or')
            %         pause(0.1)
        end
        
        FinalCell(LabelCell==CellID) = NuclID;
        
        %Associate the Nuclear tracking label Label to the corresponding Nucl
        OldNuclID = max(unique(NuclObj(NuclStruct.PixelIdxList{N})));
        NuclLabelOUT(N) = NuclLabel(OldNuclID);
        
        %Get the centroid for Cells and nucl
        NuclCentr1(N) = NuclProps(N).Centroid(1);
        NuclCentr2(N) = NuclProps(N).Centroid(2);
        
        CellCentr1(N) = CellProps(CellID).Centroid(1);
        CellCentr2(N) = CellProps(CellID).Centroid(2);
    end
    
    
end

if Debug
    SplitPerim = bwperim(CombCell);
    FocusPerim = InFocus;
    FocusPerim(SplitPerim==1) = 1;
    figure(200)
    imagesc(FinalCell)
end


%% create cytoplams object by excuding nucleus from cell obejct


%transfer Nucl image to NuclGrow
NuclGrow = LabelNucl;
%Set nuclei to 1
NuclGrow(NuclGrow>0) = 1;
%Dilate image by NuclCytoPix distance
SE = strel('disk',NuclCytoPix);
NuclGrow = imdilate(NuclGrow,SE);

%Set final Cyto image to cell image
FinalCyto = FinalCell;
%Set enlarged nuclei pixels to 0
FinalCyto(NuclGrow>0) = 0;


if Debug
    
    figure(205)
    imagesc(FinalCyto)
end

%% Save Image and measurements to VAR


AroundCellOut = Var.Analysis.AroundCellOut{CallNum};
AroundNuclOut = Var.Analysis.AroundNuclOut{CallNum};
AroundCytoOut = Var.Analysis.AroundCytoOut{CallNum};

Var.Measurements.(AroundNuclOut).SegLabel = NuclLabelOUT;
Var.Measurements.(AroundCellOut).SegLabel = NuclLabelOUT;
Var.Measurements.(AroundCytoOut).SegLabel = NuclLabelOUT;

%Transfer Centroid to Var
Var.Measurements.(AroundNuclOut).CenterX = NuclCentr1;
Var.Measurements.(AroundNuclOut).CenterY= NuclCentr2;

Var.Measurements.(AroundCellOut).CenterX= CellCentr1;
Var.Measurements.(AroundCellOut).CenterY= CellCentr2;

Var.Measurements.(AroundCytoOut).CenterX= CellCentr1;
Var.Measurements.(AroundCytoOut).CenterY= CellCentr2;


% Update objects images
Var.Img.(AroundNuclOut) = LabelNucl;
Var.Img.(AroundCellOut) = FinalCell;
Var.Img.(AroundCytoOut) = FinalCyto;




%% Diplay figure

if strcmp(Var.Figure.Display, 'on')
    
    FigNum = find(strcmp(Var.Figure.List, 'SegPhi'));
    figure(FigNum(CallNum))
    
    % RGB image fromintensity image
    BlueIntImg = (NuclImg-min(NuclImg(:)))./(max(NuclImg(:))-min(NuclImg(:)));
    BFIntImg = (InFocus-min(InFocus(:)))./(max(InFocus(:))-min(InFocus(:)));
    FluoImg(:,:,3) = BFIntImg + BlueIntImg;
    FluoImg(:,:,2) = BFIntImg;
    FluoImg(:,:,1) = BFIntImg;
    FluoImg(FluoImg>1) = 1;
    subplot(2,2,1); image(FluoImg); title('Input Images');
    
    
    P_Cell = bwperim(FinalCell);
    P_Nucl = bwperim(LabelNucl);
    P_OFF = bwperim(CombOFFCell);
    PerimImg(:,:,1) = BFIntImg+P_Cell;
    PerimImg(:,:,2) = BFIntImg+P_Nucl;
    PerimImg(:,:,3) = BFIntImg+P_OFF;
    PerimImg(PerimImg>1) = 1;
    subplot(2,2,2); image(PerimImg); title('Nucl and Cell on BF');
    subplot(2,2,3); imagesc(FinalCell); title([num2str(NumCell), ' Cells segmented from ', num2str(StartNucl), ' initial Nuclei']);
    %Generate segmentation image with cyto in green and nulceus in blue
    CombinedImg(:,:) = zeros(size(LabelNucl));
    CombinedImg(:,:,2) = im2bw(FinalCell)-im2bw(double(LabelNucl));
    CombinedImg(:,:,3) = im2bw(double(LabelNucl));
    
    subplot(2,2,4); image(CombinedImg); title('Finals Cell around Nucl Image');
    
    if Debug
        figure(803)
        image(PerimImg); title('Nucl and Cell on BF');
        FileName = [Var.Analysis.OutPath,'_NuclBFSeg.png'];
        %  imwrite(PerimImg,FileName,'png')
    end
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;





%%

function ImOut = CutWithBorder(ImIn, Corners, Border)
%corner: [Ymin, Ymax, Xmin, Xmax]
ImMid = ImIn(Corners(1):Corners(2),Corners(3):Corners(4));

ImOut = zeros(size(ImMid,1)+2*Border,size(ImMid,2)+2*Border);

ImOut(Border+1:size(ImMid,1)+Border,Border+1:size(ImMid,2)+Border)  = ImMid;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  REGION SEG                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seg = region_seg(I,init_mask,max_its,alpha,display)
% Region Based Active Contour Segmentation
%
% seg = region_seg(I,init_mask,max_its,alpha,display)
%
% Inputs: I           2D image
%         init_mask   Initialization (1 = foreground, 0 = bg)
%         max_its     Number of iterations to run segmentation for
%         alpha       (optional) Weight of smoothing term
%                       higer = smoother.  default = 0.2
%         display     (optional) displays intermediate outputs
%                       default = true
%
% Outputs: seg        Final segmentation mask (1=fg, 0=bg)
%
% Description: This code implements the paper: "Active Contours Without
% Edges" By Chan Vese. This is a nice way to segment images whose
% foregrounds and backgrounds are statistically different and homogeneous.
%
% Example:
% img = imread('tire.tif');
% m = zeros(size(img));
% m(33:33+117,44:44+128) = 1;
% seg = region_seg(img,m,500);
%
% Coded by: Shawn Lankton (www.shawnlankton.com)
%------------------------------------------------------------------------

%-- default value for parameter alpha is .1
if(~exist('alpha','var'))
    alpha = .2;
end
%-- default behavior is to display intermediate outputs
if(~exist('display','var'))
    display = true;
end
%-- ensures image is 2D double matrix
I = im2graydouble(I);

%-- Create a signed distance map (SDF) from mask
phi = mask2phi(init_mask);
FillSeg = zeros(size(I));
%--main loop
for its = 1:max_its   % Note: no automatic convergence test
    idx = find(phi <= 1.2 & phi >= -1.2);  %get the curve's narrow band
    
    
    
    %-- find interior and exterior mean
    upts = find(phi<=0);             % interior points
    vpts = find(phi>0);                  % exterior points
    u = sum(I(upts))/(length(upts)+eps); % interior mean
    v = sum(I(vpts))/(length(vpts)+eps); % exterior mean
    
    F = (I(idx)-u).^2-(I(idx)-v).^2;         % force from image information
    
    curvature = get_curvature(phi,idx);  % force from curvature penalty
    
    dphidt = F./max(abs(F)) + alpha*curvature;  % gradient descent to minimize energy
    
    %-- maintain the CFL condition
    %     Sort = sort(dphidt,'descend');
    %     if length(Sort)<10
    %         Max10 = mean(Sort);
    %     else
    %         Max10 = mean(Sort(1:10));
    %     end
    %     dt = .45/(Max10+eps);
    
    dt = .45/(max(dphidt)+eps);
    %%ADDED from original code
    %Prevents error if dt is too large
    MaxDT = 5;
    if abs(dt)> MaxDT
        dt = sign(dt)*MaxDT;
    end
    %-- evolve the curve
    try
        phi(idx) = phi(idx) + dt.*dphidt;
    catch
        warning('Limit evolve Stuck')
        break
        %             S_phi = size(phi(idx))
        %             S_dt = size(dt)
        %             SDp = size(dphidt)
        %             imagesc(seg)
        %             error('bla')
        
    end
    
    %-- Keep SDF smooth
    phi = sussman(phi, .5);
    
    
    
    seg = phi<=0; %-- Get mask from levelset
    
    PrevSeg = FillSeg;
    FillSeg = imfill(seg, 'holes');
    DiffP(its) = abs(sum(FillSeg(:)-PrevSeg(:)));
    if its >100 && mean(DiffP(end-10:end))== 0
        break
    end
    
    %     %-- intermediate output
    %     if((display>0)&&(mod(its,20) == 0))
    %       showCurveAndPhi(I,phi,its);
    %     end
end

%-- final output
%   if(display)
%     showCurveAndPhi(I,phi,its);
%   end

%-- make mask from SDF
seg = phi<=0; %-- Get mask from levelset


%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------


%-- Displays the image with curve superimposed
function showCurveAndPhi(I, phi, i)
imshow(I,'initialmagnification',200,'displayrange',[0 255]); hold on;
contour(phi, [0 0], 'g','LineWidth',4);
contour(phi, [0 0], 'k','LineWidth',2);
hold off; title([num2str(i) ' Iterations']); drawnow;

%-- converts a mask to a SDF
function phi = mask2phi(init_a)
phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;

%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
[dimy, dimx] = size(phi);
[y x] = ind2sub([dimy,dimx],idx);  % get subscripts

%-- get subscripts of neighbors
ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

%-- bounds checking
ym1(ym1<1) = 1; xm1(xm1<1) = 1;
yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;

%-- get indexes for 8 neighbors
idup = sub2ind(size(phi),yp1,x);
iddn = sub2ind(size(phi),ym1,x);
idlt = sub2ind(size(phi),y,xm1);
idrt = sub2ind(size(phi),y,xp1);
idul = sub2ind(size(phi),yp1,xm1);
idur = sub2ind(size(phi),yp1,xp1);
iddl = sub2ind(size(phi),ym1,xm1);
iddr = sub2ind(size(phi),ym1,xp1);

%-- get central derivatives of SDF at x,y
phi_x  = -phi(idlt)+phi(idrt);
phi_y  = -phi(iddn)+phi(idup);
phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
    +0.25*phi(iddr)+0.25*phi(idul);
phi_x2 = phi_x.^2;
phi_y2 = phi_y.^2;

%-- compute curvature (Kappa)
curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
    (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);

%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)
[dimy, dimx, c] = size(img);
if(isfloat(img)) % image is a double
    if(c==3)
        img = rgb2gray(uint8(img));
    end
else           % image is a int
    if(c==3)
        img = rgb2gray(img);
    end
    img = double(img);
end

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
% forward/backward differences
a = D - shiftR(D); % backward
b = shiftL(D) - D; % forward
c = D - shiftD(D); % backward
d = shiftU(D) - D; % forward

a_p = a;  a_n = a; % a+ and a-
b_p = b;  b_n = b;
c_p = c;  c_n = c;
d_p = d;  d_n = d;

a_p(a < 0) = 0;
a_n(a > 0) = 0;
b_p(b < 0) = 0;
b_n(b > 0) = 0;
c_p(c < 0) = 0;
c_n(c > 0) = 0;
d_p(d < 0) = 0;
d_n(d > 0) = 0;

dD = zeros(size(D));
D_neg_ind = find(D < 0);
D_pos_ind = find(D > 0);
dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
    + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
    + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

D = D - dt .* sussman_sign(D) .* dD;

%-- whole matrix derivatives
function shift = shiftD(M)
shift = shiftR(M')';

function shift = shiftL(M)
shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
shift = shiftL(M')';

function S = sussman_sign(D)
S = D ./ sqrt(D.^2 + 1);


