function Var = TrackObjectsSingle(Var, CallNum)
tic
Debug = 0;      %Set Debug to 0 to prevent image display

if nargin == 1
    CallNum = 1;
end

%% Get Track Object and Tracking distance from Var
TrackObj = Var.Analysis.TrackObj{CallNum};
MaxDistance = Var.Analysis.TrackMaxDistance{CallNum};
TrackBack = 5; %Number of frames to compare in the past

%Read data from 1st frame
F = 1;
%Read Frame.mat file
% ImgFolder = fileparts(Var.Analysis.OutPath);
% FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(F,'%05d'),'.mat']);
% load(FileName)

%Get X and Y centers from measurement
PresentCenter = [Var.SegMeas(F).(TrackObj).CenterX, Var.SegMeas(F).(TrackObj).CenterY];

TrackLabel{F} = Var.SegMeas(F).(TrackObj).SegLabel;
MaxLabel = max(Var.SegMeas(F).(TrackObj).SegLabel);



for F = 2:Var.Analysis.NumFrame
    %Set coordinates from present center to previous center
    PreviousCenter = PresentCenter;
    
    %Get X and Y centers from measurement
    PresentCenter = [Var.SegMeas(F).(TrackObj).CenterX, Var.SegMeas(F).(TrackObj).CenterY];
    
    if ~isempty(PreviousCenter)
        
        
        if Debug ==1
            assignin('base','PresentCenter',PresentCenter);
            assignin('base','PreviousCenter',PreviousCenter);
        end
        
        % runs a simple routine to find the distance between the reference
        % centers and the image center
        %DistanceMap is an array(NumberReferenceXNumberCenter) of all the
        %distance calculated, CellPos is the Number of the object center closest to the reference center(given by the indice)
        %MinimalDistance is the smallest distance calculated between the
        %reference center and the Image center
        [MinimalDistance, CellPos, DistanceMap] = findclosest(PresentCenter, PreviousCenter);
        
        %Set CellPos to zero for minimal distance larger then TrackDistance
        CellPos(MinimalDistance>MaxDistance) = 0;
        
        if Debug ==1
            assignin('base','DistanceMap',DistanceMap);
            assignin('base','MinimalDistance',MinimalDistance);
            assignin('base','Orig_CellPos',CellPos);
        end
        
        NumberCenter = size(PresentCenter,1);
        PrevNumberCenter = size(PreviousCenter,1);
        for m = 1:PrevNumberCenter
            
            % fprintf(num2str(m))
            CellNum = find(CellPos==m);
            
            
            %if number of CellPos equal to an object is larger than one =>
            %there are duplicate matches
            if length(CellNum) > 1
                %sort the duplicate match by distance order
                [SortMD, IndiceSortMD] = sort(MinimalDistance(CellNum), 'ascend');
                %The center with the closest distance is kept
                %Tries to find if the next best match for the other centers
                for k = 2:length(CellNum)
                     %find all the other Object numbers not equal to the duplicate match
                    OppCellNum = find(CellPos(:)~=m);
                    %sort distance match to find next best center position
                    [SortMap, IndiceSortMap] = sort(DistanceMap(CellNum(IndiceSortMD(k)),:));
                    
                    DefinedLabel = 0;
                    iter = 1;
                    while DefinedLabel == 0
                        %Starts iter from 2 because we know that the closest object
                        %is already assigned
                        iter = iter+1;
                        % Check if distance is smaller than maximal distance
                        % movement allowed
                        if iter < length(SortMap) && SortMap(iter) < MaxDistance
                            %Test if the center is already in the list of
                            %matched centers
                            if isempty(find(CellPos(OppCellNum) == IndiceSortMap(iter)))
                                %updates CellPos and MinimalDistance with the new
                                %values
                                CellPos(CellNum(IndiceSortMD(k))) = IndiceSortMap(iter);
                                MinimalDistance(CellNum(IndiceSortMD(k))) = SortMap(iter);
                                DefinedLabel = 1;
                            end
                        else
                            %No Position found Set CellPos to 0 and
                            %MinDistance to Inf
                            CellPos(CellNum(IndiceSortMD(k))) = 0;
                            MinimalDistance(CellNum(IndiceSortMD(k))) = Inf;
                            DefinedLabel = 1;
                            
                        end
                    end
                end
            end
        end % end loop through center search
        
        
        
        if Debug ==1
            assignin('base','Final_CellPos',CellPos);
        end
    else %WHAT happens if no previous center???
        % Set CellPos to 0
        NumberCenter = size(PresentCenter,1);
        CellPos = zeros(NumberCenter,1);
        
        % and mindist to Inf
        MinimalDistance = ones(NumberCenter,1).*Inf;
    end
    
    
    %% Compare with previous time points
    
    %     if isfield(Var.Measurements.(Var.Analysis.TrackObj{CallNum}), 'Label') &&...
    %             length(Var.Measurements.(Var.Analysis.TrackObj{CallNum}).Label) >= Var.Analysis.FrameIter-1
    %         %Load previous labels
    %         PreviousLabel = Var.Measurements.(Var.Analysis.TrackObj{CallNum}).Label{Var.Analysis.FrameIter-1};
    %     else
    %         PreviousLabel = [];
    %     end
    
    
    %find unassigned present nuclei
    NoMatchIndex = find(CellPos == 0);
    assignin('base','NoMatchIndex',NoMatchIndex);
    if ~isempty(NoMatchIndex)
        NumberNoMatch = length(NoMatchIndex);
        FoundMatch = 0;
        ComparedLabel = TrackLabel{F-1};
        %Loop through previous frames to check the past labels
        iter = 2;
        while (F-iter) > 0 && iter <= TrackBack
            % F = F
            % iter = iter
            BestMatchDist = [];
            BestMatchLabel = [];
            %PastLabel = Var.Measurements.(Var.Analysis.TrackObj{CallNum}).Label{Var.Analysis.FrameIter-iter};
            
            %Load previous labels
            PastLabel = TrackLabel{F-iter};
            
            
            %find label present in the past but not in Frame-1 matrix
            [MissingLabels, MissIndex] = setdiff(PastLabel, ComparedLabel);
            %             MissingLabels = MissingLabels
            %                      MissIndex = MissIndex
            
            
            %Check whether previous cells in other frames were not assigned
            if ~isempty(MissingLabels) && ~isempty(NoMatchIndex)
                %Get X and Y centers from measurement
                MissCenter = [Var.SegMeas(F-iter).(TrackObj).CenterX(MissIndex), Var.SegMeas(F-iter).(TrackObj).CenterY(MissIndex)];
                %Get centers from unassigned past labels
                %Get centers from present unmatched cells
                NoMatchCenter = [PresentCenter(NoMatchIndex,1),PresentCenter(NoMatchIndex,2)];
                %Calculate minimal distance array
                [MinimalDistance, CellPos_bis, DistanceMap] = findclosest(NoMatchCenter, MissCenter);
                %Check if the distance is small enough if distance too large
                %set cellpos to 0
                %MinimalDistance = MinimalDistance
                CellPos_bis(MinimalDistance > MaxDistance/1) = 0;      %MaxDistance/4???
                
                
                
                %Loop through unmatched items
                for m = 1:NumberNoMatch
                    if CellPos_bis(m) ~= 0
                        UniqueLabel = find(CellPos_bis == CellPos_bis(m));
                        if length(UniqueLabel) >1
                            
                            [~, SortUnique] = sort(MinimalDistance(UniqueLabel), 'ascend');
                            CellPos_bis(UniqueLabel(SortUnique(2:end))) = 0;
                            MinimalDistance(UniqueLabel(SortUnique(2:end))) = +inf;
                            
                        end
                        
                    end
                    if CellPos_bis(m) ~= 0
                        FoundMatch = 1;
                        BestMatchDist(iter,m) = MinimalDistance(m);
                        BestMatchLabel(iter,m) = MissingLabels(CellPos_bis(m));
                    end
                end
            end
            ComparedLabel = union(ComparedLabel, MissingLabels);
            iter = iter +1;
        end
        if FoundMatch == 1
            %Find in history shorter distance
            BestMatchDist(BestMatchDist == 0) = +inf;
            %     BestMatchDist = BestMatchDist
            %     BestMatchLabel = BestMatchLabel
            for m = 1:size(BestMatchDist,2)
                [SortBD, IndiceSortBD] = sort(BestMatchDist(:,m));
                BestLabel = BestMatchLabel(IndiceSortBD(1),m);
                CellPos(NoMatchIndex(m)) = -1*BestMatchLabel(IndiceSortBD(1),m);
            end
        end
    end
    
    
    % assignin('base','PreviousLabel',PreviousLabel);
    % initalize matrices
    
    TrackLabel{F} = zeros(NumberCenter,1);
    for m = 1:NumberCenter
        if CellPos(m) > 0
            TrackLabel{F}(m) = TrackLabel{F-1}(CellPos(m));
        elseif CellPos(m) < 0
            TrackLabel{F}(m) = -1*(CellPos(m));
        elseif CellPos(m) == 0
            MaxLabel = MaxLabel +1;
            TrackLabel{F}(m) = MaxLabel;
        end
        
    end
    assignin('base','TrackLabel',TrackLabel);
    
    
    
    %     CheckProps = regionprops(PresentImg, 'Area', 'PixelIdxList', 'Centroid', 'EquivDiameter');
    %     CheckImg = zeros(size(NewImg));
    %     NumObjects = size(CheckProps,1);
    %     ObjName = 'Nucl';
    %     iter = 0;
    %    for i = 1:NumObjects
    %        CheckImg = zeros(size(NewImg));
    %         if CheckProps(i).Area == 0
    %              fprintf(['No ', ObjName, ' # ', num2str(i)  , '\n'])
    %             %object absent from Frame
    %         else
    %             iter = iter +1;
    %             %Object Present in Frame
    %             %Check if label and Object Number are not equivalent
    %             if Label(iter) ~=  NewImg(round(CheckProps(i).Centroid(2)), round(CheckProps(i).Centroid(1)));
    %                 fprintf(['Wrong Label for ', ObjName, ' # ', num2str(i)  , ' with Label ', num2str(Label(iter)), '\n'])
    %             else
    %                  fprintf(['OKay Label for ', ObjName, ' # ', num2str(i)  , ' with Label ', num2str(Label(iter)), '\n'])
    %             end
    %             CheckImg(CheckProps(i).PixelIdxList) = i;
    %             figure(104)
    %             imagesc(CheckImg)
    %             pause(0.5)
    %         end
    %    end
    
    
    %% Display Image
    
    %     %Display Image from combined objects
    %     if strcmp(Var.Figure.Display, 'on')
    %         FigNum = find((strcmp(Var.Figure.List, 'TrackObj')));
    %         figure(FigNum(CallNum))
    %         plot(PresentCenter(:,1),PresentCenter(:,2), '.b')
    %         hold on
    %         plot(PreviousCenter(:,1),PreviousCenter(:,2), '.r')
    %         axis([1 Var.Analysis.ImgSize(2) 1 Var.Analysis.ImgSize(1)])
    %         set(gca, 'YDir', 'reverse')
    %
    %         for m = 1:size(PresentCenter,1)
    %             text(PresentCenter(m,1),PresentCenter(m,2), num2str(TrackLabel{F}(m)), 'FontSize', 16);
    %         end
    %         for m = 1:size(PreviousCenter,1)
    %             text(PreviousCenter(m,1),PreviousCenter(m,2), num2str(TrackLabel{F-1}(m)),'FontAngle', 'italic', 'FontSize', 16); %
    %         end
    %         hold off
    %
    %         title(['Frame: ',num2str(F), ' MaxLabel: ', num2str(MaxLabel)])
    %
    %     end
    %     drawnow
    assignin('base', 'TrackLabel', TrackLabel)
    if strcmp(Var.Figure.Display, 'on')
        FigNum = find((strcmp(Var.Figure.List, 'TrackObj')));
        
        figure(FigNum(CallNum))
        hold on
        
        
        for C = 1:length(TrackLabel{F})
            
            RGBcolor = hsv2rgb([(F-1)/Var.Analysis.NumFrame,1,1]);
            % C = C
            %             TrackLabel{F}(C)
            PrevC = find(TrackLabel{F-1} == TrackLabel{F}(C));
            if ~isempty(PrevC)
                plot([PreviousCenter(PrevC,1),PresentCenter(C,1)],[PreviousCenter(PrevC,2),PresentCenter(C,2)], ...
                    'Color' , RGBcolor,...
                    'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', RGBcolor)
            else
                plot(PresentCenter(C,1),PresentCenter(C,2), 'Color' , RGBcolor, 'Marker', 'o', 'MarkerSize', 2)
                text(PresentCenter(C,1),PresentCenter(C,2), num2str(TrackLabel{F}(C)), 'FontSize', 12);
            end
        end
        if F ==2
            for C = 1:length(TrackLabel{F-1})
                plot(PreviousCenter(C,1),PreviousCenter(C,2), 'Color' , RGBcolor, 'Marker', 'o', 'MarkerSize', 4)
                
                text(PreviousCenter(C,1),PreviousCenter(C,2), num2str(TrackLabel{F-1}(C)), 'FontSize', 12);
            end
            
        end
        axis([1 Var.Analysis.ImgSize(2) 1 Var.Analysis.ImgSize(1)])
        set(gca, 'YDir', 'reverse')
        
        %         for m = 1:size(PresentCenter,1)
        %             text(PresentCenter(m,1),PresentCenter(m,2), num2str(TrackLabel{F}(m)), 'FontSize', 16);
        %         end
        %         for m = 1:size(PreviousCenter,1)
        %             text(PreviousCenter(m,1),PreviousCenter(m,2), num2str(TrackLabel{F-1}(m)),'FontAngle', 'italic', 'FontSize', 16); %
        %         end
        hold off
        
        title(['Frame: ',num2str(F), ' MaxLabel: ', num2str(MaxLabel)])
        
    end
    drawnow
    %     pict_name = ['Track3_',num2str(Var.Analysis.FrameIter,'%02d')];
    %     print('-dpng',FigNum(CallNum),pict_name);
    % if strcmp(Var.Analysis.FrameAnalyzed,'Last')
    %     Var.Analysis.MovieObj{10} = close(Var.Analysis.MovieObj{10});
    % end
    
end


%% Save TrackLabel to Var structure

Var.TrackLabel = TrackLabel;


%% link present object to other object in the analysis
% if isfield(Var.Analysis, 'LinkedObj') && ~isempty(Var.Analysis.LinkedObj{CallNum})
%     LinkObjName = Var.Analysis.LinkedObj{CallNum};
%     LinkObjImg = Var.Img.(LinkObjName);
%     for i = 1:length(ObjProps)
%       %  i = i
%         %Get Linked object number from centroid of present object
%         %         SumImg = LinkObjImg; SumImg(ObjProps(i).PixelIdxList) =
%         %         SumImg(ObjProps(i).PixelIdxList)+10; figure(120)
%         %         imagesc(SumImg)
%         % %         hold on %
%         % plot(round(ObjProps(i).Centroid(1)),round(ObjProps(i).Centroid(2)),
%         % 'rx') % %         hold off
%         %         pause(1)
%
%         if Debug == 1
%             figure(200)
%             LOI = zeros(size(LinkObjImg));
%             LOI(LinkObjImg>0) = 1;
%             LOI(ObjProps(i).PixelIdxList) = 2;
%             imagesc(LOI)
%             pause(1)
%         end
%         LinkObjNum = unique(LinkObjImg(ObjProps(i).PixelIdxList));
%         LinkObjNum = LinkObjNum(LinkObjNum>0);
%       %  LinkObjNum = LinkObjNum
%         %%PB finding cells outside of main object Quick fix: set linked
%         %%object label to 0 Should be fixed in the segmentation routine...
%         if length(LinkObjNum) == 1;
%             %LinkObjNum =
%             %LinkObjImg(round(ObjProps(i).Centroid(2)),round(ObjProps(i).Centroid(1)))
%             LinkedLabel = Var.SegMeas(F).(LinkObjName).Label{Var.Analysis.FrameIter}(LinkObjNum);
%             Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).ParentLabel{Var.Analysis.FrameIter}(i) = LinkedLabel;
%        elseif length(LinkObjNum) > 1;
%            AllLinkObjNum = LinkObjImg(ObjProps(i).PixelIdxList);
%            AllLinkObjNum = AllLinkObjNum(AllLinkObjNum>0);
%            NumEach = zeros(length(LinkObjNum),1);
%            for LO = 1:length(LinkObjNum)
%                NumEach(LO) = length(find(AllLinkObjNum == LinkObjNum(LO)));
%            end
%            [MaxPix , IndMax] = max(NumEach);
%            LinkedLabel = Var.SegMeas(F).(LinkObjName).Label{Var.Analysis.FrameIter}(LinkObjNum(IndMax));
%             Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).ParentLabel{Var.Analysis.FrameIter}(i) = LinkedLabel;
%         elseif length(LinkObjNum) == 0;
%             Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).ParentLabel{Var.Analysis.FrameIter}(i) = 0;
%         end
%     end
% end


% if Debug == 1 && CallNum == 2
%     Hue = linspace(0,0.9,30);
%     if strcmp(Var.Analysis.FrameAnalyzed, 'First')
%         mov_name = ['/Users/serge/Projects/ImageAnalysis/Data/130614/0614_Ctt1B_0p2M/TestObjectTrack.avi'];
%         Var.aviobj = avifile(mov_name, 'fps', 2, 'Quality', 100, 'Compression', 'None');
%     end
%
%     ObjectImageH = zeros(size(Var.Img.(Var.Analysis.TrackObj{CallNum})));
%     ObjectImageS = zeros(size(Var.Img.(Var.Analysis.TrackObj{CallNum})));
%     ObjectImageV = (Var.Img.GFPImg-min(Var.Img.GFPImg(:)))./(max(Var.Img.GFPImg(:))-min(Var.Img.GFPImg(:))); %ones(size(Var.Img.(Var.Analysis.TrackObj{CallNum})));
%     for j = 1:Var.Analysis.MaxLabel(CallNum)
%         HueInd = rem(j-1,length(Hue))+1;
%         ObjectImageH(Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).PixelList{Var.Analysis.FrameIter,(Label == j)}) = Hue(HueInd);
%         ObjectImageS(Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).PixelList{Var.Analysis.FrameIter,(Label == j)}) = 1;
%         ObjectImageB(Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).PixelList{Var.Analysis.FrameIter,(Label == j)}) = 1;
%     end
%     FigureNb = 432;
%     figure(FigureNb)
%     pos = [180; 251; 896; 672];
%     set(FigureNb, 'Position', pos)
%     set(FigureNb, 'renderer', 'OpenGL')
%     HSVObjectImage(:,:,1)  = ObjectImageH;
%     HSVObjectImage(:,:,2)  = ObjectImageS;
%     HSVObjectImage(:,:,3)  = ObjectImageV; %label2rgb(ObjectImage);
%
%     RGBObjectImage = hsv2rgb(HSVObjectImage);
%     image(RGBObjectImage)
%
%     Center = [Var.SegMeas(F).(Var.Analysis.TrackObj{CallNum}).CenterX{Var.Analysis.FrameIter}, Var.Measurements.(Var.Analysis.TrackObj{CallNum}).CenterY{Var.Analysis.FrameIter}];
%     for m = 1:length(Center)
%         text(Center(m,1),Center(m,2), num2str(Var.Measurements.(Var.Analysis.TrackObj{CallNum}).Label{Var.Analysis.FrameIter}(m)), 'FontSize', 14,'Color', 'r'); %
%     end
%
%     set(gca,'ytick',[], 'xtick', []);
%     pause(0.1)
%    Frame = getframe(FigureNb);
%   Var.aviobj = addframe(Var.aviobj,Frame);
%
%     if strcmp(Var.Analysis.FrameAnalyzed, 'Last')
%         Var.aviobj = close(Var.aviobj);
%     end
%
% end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;


%%
%%%%%%%%%%%%%%%%%%%%%%%
%      SUB function   %
%%%%%%%%%%%%%%%%%%%%%%%


function [MinDist,CellPos, DistanceRefCenter] = findclosest(ReferenceCenter, ImageCenter)

NumberReference = size(ReferenceCenter,1);
NumberObjectImage = size(ImageCenter,1);
DistanceRefCenter = zeros(NumberReference, NumberObjectImage);

for i = 1:NumberReference
    for j = 1:NumberObjectImage
        DistanceRefCenter(i,j) = sqrt((ReferenceCenter(i,1)-ImageCenter(j,1))^2 ...
            + (ReferenceCenter(i,2)-ImageCenter(j,2))^2);
    end
end
[MinDist,CellPos] = min(DistanceRefCenter,[],2);