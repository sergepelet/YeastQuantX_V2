function Thumb = DisplayCellImage(Path, Position, CellLabel)
%% Display image from cells based on their label and position from the YeastQuant analysis
% Path: can either be the access path to the Data.mat file : /Volumes/DMF/GROUPS/gr_Pelet/Serge/2014/141211/1211_ySP569ButShort/img_Pos01_Data.mat
%           or the folder where the Data.mat file is located     /Volumes/DMF/GROUPS/gr_Pelet/Serge/2014/141211/1211_ySP569ButShort/
            % !! should end with /
% Positio: Number from position in matlab analysis

% CellLabel: Label number from cell to be displayed


%% Set Time point and illumniation to display
TimePoint =  [1,8,14]; % [1,7,16]; %Mpk1Zymp [1,8,15];  %Ste7 %[1,8,15]; %Ste7_CC [1,7,16]; %Mpk1Zymp      %[1,8,14]; %Ste7_yox1

% Provide string to identify illumination to be loaded. It has to match the
% filter settings specified in yeastquant datadase
% Set Color for image display : R, G, B or a combination of the letters.
% Illum = {'YFP', 'CFP', 'RFP'};
% Color = {'G', 'B', 'R'};
% 
% Illum = {'BF0', 'CFP', 'RFP'};
% Color = {'RGB', 'B', 'R'};
% 
% Illum = {'BF0', 'YFP', 'RFP'};
% Color = {'RGB', 'G', 'R'};

Illum = {'BF0','CFP', 'RFP', 'YFP'};
Color = {'RGB','B', 'R', 'G'};

%Reference data to find Object in Data.mat file
RefObject = 'Nucl';
RefIllum = 'CorrCFP';
CheckLabel = CellLabel + Position/1000;

%Set object displayed on image and on wich illumination
%Leave empty '' if no object has to be displayed
%At least one object has to be specified to retrieve the coordinates of the
%object
Object = {'Cell', '', '',  ''};
%Define color of object: set to 1 for white and 0 for black 
ObjectColor = {0,[],[], []}; 

%Set intensity limits for each image type
%IntLimit = {[50,300],[50,1500],[300,2500] };
%IntLimit = {[500,6000],[50,1000],[300,3000] };  %Ste7DS ??
%IntLimit = {[5000,20000],[500,2000],[300,2000] };  %Mpk1Sensor
% IntLimit = {[2200,8500],[100,1300],[100,800] };  %Mpk1Sensor

IntLimit = {[5000,20000],[100 1500],[300,4000], [100,500]  }; %Whi5 Ste7DS

%Size of thumbnail image
ThumbSize = 60;

%% File loading
%Get Folder name and extension
[Folder, ~, Extension] = fileparts(Path);
if strcmp(Extension, '.mat')
    %Path is the direct access path to the Data.mat file
    VarFilePath = Path;

else
    VarFilePath = fullfile(Path, ['Pos', num2str(Position, '%02d_Data.mat')]);
end
%Load Var from Data.mat file
load(VarFilePath)

%Loop through all time points
for T = 1:length(TimePoint)
    %Get index of cell in measurement label matrix
    CellInd = find(Var.Measurements.(RefObject).(RefIllum).CheckLabel(TimePoint(T),:) == CheckLabel);
    %Get XY center of object
    Xpos = Var.Measurements.(RefObject).(RefIllum).CenterX(TimePoint(T),CellInd);
    Ypos = Var.Measurements.(RefObject).(RefIllum).CenterY(TimePoint(T),CellInd);
    
    %Loop through each illumination
    for I = 1:length(Illum)
        %Match Illum string to the Filter settings in
        IllumNum = find(strcmpi(Var.Experiment.Filter, Illum{I}));
        %Get Image file path from illumination number and time point
        RawFile = Var.Analysis.FilePath{IllumNum,TimePoint(T)};
        %Correct Image for server access path
        [RawFolder, ImgName] = fileparts(RawFile);
        [~, PosName] = fileparts(RawFolder);
        %New Image path
        FullImgPath = fullfile(Folder, PosName, [ImgName,'.tif']);
        %Read image
        Img = imread(FullImgPath);
        %Normalize image
        Img = double(Img);
        %     MinI = min(Img(:))
        %     MaxI = max(Img(:))
        if isempty(IntLimit{I})
            Img = (Img-min(Img(:)))./(max(Img(:))-min(Img(:)));
        else
            Img = (Img-IntLimit{I}(1))./(IntLimit{I}(2)-IntLimit{I}(1));
            Img(Img>1) = 1;
            Img(Img<0) = 0;
        end
        
        %Get pixels from object
        if ~isempty(Object{I})
            PixelList = Var.Measurements.(Object{I}).PixelList{TimePoint(T),CellLabel};
        else
            PixelList = [];
        end
        %Create Object image
        ObjImg = zeros(size(Img));
        ObjImg(PixelList) =1;
        %%get perimeter of object
        
        %Create thumbnail image and object
        %Calculate X Y coordinates of thumbnail based on position of the
        %image
        Crop = [max(round(Ypos-ThumbSize/2),1), min(round(Ypos+ThumbSize/2)-1,size(Img,1)), ...
            max(round(Xpos-ThumbSize/2),1), min(round(Xpos+ThumbSize/2)-1,size(Img,2))];
        
        %Crop image
        Thumb = Img(Crop(1):Crop(2),Crop(3):Crop(4));
        %Crop Object image
        ThumObj = ObjImg(Crop(1):Crop(2),Crop(3):Crop(4));
        
        %Get perimeter of object
        ThumbPerim = bwperim(ThumObj);
        if ~isempty(Object{I})
            %Set pixels belonging to perimiter of object in balck or white
            Thumb(ThumbPerim == 1) = ObjectColor{I};
        end
        
        %Create RGB image and set each color to Thumb image if needed
        RGBImg = zeros(size(Thumb,1), size(Thumb,2),3);
        if ~isempty(strfind(Color{I}, 'R'))
            RGBImg(:,:,1) = Thumb;
        end
        if ~isempty(strfind(Color{I}, 'G'))
            RGBImg(:,:,2) = Thumb;
        end
        if ~isempty(strfind(Color{I}, 'B'))
            RGBImg(:,:,3) = Thumb;
        end

        
        %Create empty image containing all time points and illumination
        if I == 1
            if T == 1
                TimeLapseImg = zeros(size(Thumb,1)*length(Illum),size(Thumb,2)*length(TimePoint),3);
            end
        end

        %Set sub-matrix of TimeLapseImg to RGB
        TimeLapseImg((I-1)*size(Thumb,1)+1:I*size(Thumb,1),(T-1)*size(Thumb,2)+1:T*size(Thumb,2),:) = RGBImg;
        
    end
end

%Display image
figure(Position*1000+CellLabel)
imshow(TimeLapseImg)
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [length(Illum)*6 length(TimePoint)*6]);
set(gcf, 'PaperPosition', [0 0 length(Illum)*6 length(TimePoint)*6]);



