function Thumb = CheckLabel2TIFF_Zplanes(Folder, CheckLabelList)

%Size of thumbnail image
ThumbSize = 60;
%Reference data to find Object in Data.mat file
RefObject = 'Cell';
RefIllum = 'CorrGFP';



for C = 1:length(CheckLabelList)
    
    Position(C) = round(rem(CheckLabelList(C),1)*1000);
    CellLabel(C) = floor(CheckLabelList(C));
    
    ThumbFolder{C} = fullfile(Folder, ['Thumb_', num2str(Position(C),'%03d_'), num2str(CellLabel(C),'%05d')]);
    mkdir(ThumbFolder{C})
end

PosList = unique(Position);
PosIter = 0;
for P = PosList
    PosIter = PosIter+1;
    %% Data.mat File loading
    VarFilePath = fullfile(Folder, 'DATAOut', ['Pos_', num2str(PosList(PosIter), '%03d_Data.mat')]);
    
    %Load Var from Data.mat file
    load(VarFilePath)
    
    NumTpts = size(Var.Analysis.FilePath,2);
    NumIllum = size(Var.Analysis.FilePath,1);
    %% Loop through all time points
    for T = 1:NumTpts
        
        %% Loop through all Illuminations
        for I = 1:NumIllum
            
            %Get Image file path from illumination number and time point
            RawFile = Var.Analysis.FilePath{I,T};
            SkipBack = 0;
            while isempty(RawFile)
                SkipBack = SkipBack+1;
                RawFile = Var.Analysis.FilePath{I,T-SkipBack};
            end
            if iscell(RawFile)
                for Z = 1:length(RawFile)
                    
                    %Correct Image for server access path
                    [RawFolder, ImgName] = fileparts(RawFile{Z});
                    [~, PosName] = fileparts(RawFolder);
                    %New Image path
                    FullImgPath = fullfile(Folder, PosName, [ImgName,'.tif']);
                    %Read image
                    Img = imread(FullImgPath);
%                     %MaxIntProjection
%                     if Z == 1
%                         MaxImg = Img;
%                     else
%                         MaxImg = max(MaxImg, Img);
%                     end
%                     
%                 end
%                 Img = MaxImg;

            %Find all Labels from current position
            LabelInd = find(Position == PosList(PosIter));
            % Loop through all labels
            for L = 1:length(LabelInd)
                
                %Get index of cell in measurement label matrix
                CellInd = find(Var.Measurements.(RefObject).(RefIllum).CheckLabel(T,:) == CheckLabelList(LabelInd(L)));
                %Get XY center of object
                Xpos = Var.Measurements.(RefObject).(RefIllum).CenterX(T,CellInd);
                Ypos = Var.Measurements.(RefObject).(RefIllum).CenterY(T,CellInd);
                
                %Calculate X Y coordinates of thumbnail based on position of the
                %image
                Crop = [max(round(Ypos-ThumbSize/2),1), min(round(Ypos+ThumbSize/2)-1,Var.Analysis.ImgSize(1)), ...
                    max(round(Xpos-ThumbSize/2),1), min(round(Xpos+ThumbSize/2)-1,Var.Analysis.ImgSize(2))];
                %Crop the image
                Thumb = Img(Crop(1):Crop(2),Crop(3):Crop(4));
                %Generate file name
                ThumbFile = fullfile(ThumbFolder{LabelInd(L)}, [ImgName,'.tif']);
                %Save tif image
                imwrite(Thumb, ThumbFile, 'tif')
            end
                end
            else
                %Correct Image for server access path
                [RawFolder, ImgName] = fileparts(RawFile);
                [~, PosName] = fileparts(RawFolder);
                %New Image path
                FullImgPath = fullfile(Folder, PosName, [ImgName,'.tif']);
                %Read image
                Img = imread(FullImgPath);
            
            
            %Find all Labels from current position
            LabelInd = find(Position == PosList(PosIter));
            % Loop through all labels
            for L = 1:length(LabelInd)
                
                %Get index of cell in measurement label matrix
                CellInd = find(Var.Measurements.(RefObject).(RefIllum).CheckLabel(T,:) == CheckLabelList(LabelInd(L)));
                %Get XY center of object
                Xpos = Var.Measurements.(RefObject).(RefIllum).CenterX(T,CellInd);
                Ypos = Var.Measurements.(RefObject).(RefIllum).CenterY(T,CellInd);
                
                %Calculate X Y coordinates of thumbnail based on position of the
                %image
                Crop = [max(round(Ypos-ThumbSize/2),1), min(round(Ypos+ThumbSize/2)-1,Var.Analysis.ImgSize(1)), ...
                    max(round(Xpos-ThumbSize/2),1), min(round(Xpos+ThumbSize/2)-1,Var.Analysis.ImgSize(2))];
                %Crop the image
                Thumb = Img(Crop(1):Crop(2),Crop(3):Crop(4));
                %Generate file name
                ThumbFile = fullfile(ThumbFolder{LabelInd(L)}, [ImgName,'.tif']);
                %Save tif image
                imwrite(Thumb, ThumbFile, 'tif')
            end
            end
        end
    end
    
end