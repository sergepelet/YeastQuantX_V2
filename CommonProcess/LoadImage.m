function Var = LoadImage(Var, CallNum)
tic
if nargin == 1
    CallNum = 1;
end
%Set variable
IllumPara = strmatch( Var.Analysis.LoadIllum{CallNum}, Var.Experiment.Filter);
%If no match found between loadIllum and filter strings it means that a
%Z-plane information is added to loadillum string
if isempty(IllumPara)
    %Find underscore
    US = strfind(Var.Analysis.LoadIllum{CallNum}, '_');
    %Get Illum string and number
    IllumStr = Var.Analysis.LoadIllum{CallNum}(1:US-1);
    IllumPara = strmatch( IllumStr, Var.Experiment.Filter);
    %Get Z-Plane number
    NumPlane = str2double(Var.Analysis.LoadIllum{CallNum}(US+1:end));
end
Pos = Var.Analysis.CurrentPos;
Frame = Var.Analysis.CurrentFrame;
NumIllum = Var.Analysis.NumIllum;
NumPos = Var.Analysis.NumPos;

%Check if Skip frame iteration is one if yes load the image
if Var.Analysis.SkipFrame(2,IllumPara) == 1
    %set NewImage to 1
    Var.Analysis.NewImage(IllumPara) = 1;

%% load from metamorph stack with SingleStack acquisition
if strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    FrameInStack = (Frame-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
    [TIFFImg] = tiffread_serge(Var.Analysis.FilePath{1}, FrameInStack);
    [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{1});
   % if Frame == 1 && CallNum == 1
  %      fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
        FrameZeroNum = (Var.Experiment.TimeZero-1)*(NumPos*NumIllum) + Pos;
        ZeroFrameImg = tiffread_serge(Var.Analysis.FilePath{1}, FrameZeroNum );
        Var.Analysis.TimeZero = ZeroFrameImg.MM_stack(4);
  %  end
    Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = TIFFImg.MM_stack(4);
    
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(TIFFImg.data);
    %% load from metamorph stack with IterSave acquisition
elseif strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'IterSave')
    
    FrameInStack = (Pos-1)*NumIllum + IllumPara;
    [TIFFImg] = tiffread_serge(Var.Analysis.FilePath{Var.Analysis.CurrentFrame}, FrameInStack);
    [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{Var.Analysis.CurrentFrame});
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
        FrameZeroNum = (Pos-1)*NumIllum + 1;
        ZeroFrameImg = tiffread_serge(Var.Analysis.FilePath{Var.Experiment.TimeZero}, FrameZeroNum );
        Var.Analysis.TimeZero = ZeroFrameImg.MM_stack(4);
    end
    Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = TIFFImg.MM_stack(4);
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(TIFFImg.data);
elseif  strcmp(Var.Experiment.Software, 'NikonElement') && strcmp(Var.Experiment.AcquisitionType, 'TiffSerie')
    Img = imread(Var.Analysis.FilePath{Frame}, 'TIFF');
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(Img);
elseif  strcmp(Var.Experiment.Software, 'NikonElement') && strcmp(Var.Experiment.AcquisitionType, 'TiffExport')
    Img = imread(Var.Analysis.FilePath{IllumPara,Frame});
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(Img);
    [Var.Analysis.CurrentFolder, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{IllumPara, Frame});
    % try
    %  Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = GetTimeStamp(Var);
    %     catch
    %       FileProps = dir(Var.Analysis.FilePath{IllumPara, Frame});
    %       TimeVec =   datevec(FileProps(1).datenum);
    %       Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = TimeVec(4)*3600 +TimeVec(5)*60+TimeVec(6) ;
    %     end
    
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
        % try
        %       Var.Analysis.TimeZero = GetTimeStamp(Var, Var.Experiment.TimeZero);
        %         catch
        %         FileProps = dir(Var.Analysis.FilePath{IllumPara, Var.Experiment.TimeZero});
        %         TimeVecZero =   datevec(FileProps(1).datenum);
        %         Var.Analysis.TimeZero = TimeVecZero(4)*3600 +TimeVecZero(5)*60+TimeVecZero(6) ;
        %         end
    end
    
    %% load from metamorph Multidimensional data
elseif  strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
    if exist('NumPlane')
        Img = imread(Var.Analysis.FilePath{IllumPara,Frame},NumPlane);
    else
        Img = imread(Var.Analysis.FilePath{IllumPara,Frame});
    end
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(Img);
    [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{IllumPara, Frame});
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
    end
    
    %% load from MicroManager from MultiDimensionnal data
elseif strcmp(Var.Experiment.Software, 'MicroManager') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
    
    if isfield(Var.Analysis, 'Zstack') && Var.Analysis.Zstack(IllumPara) >1
        %Load Z-Stack and calculate MaxProjection
        Img = LoadZstack(Var.Analysis.FilePath{IllumPara,Frame});
        [Var.Analysis.CurrentFolder, Var.Analysis.CurrentFileName]  = fileparts(Var.Analysis.FilePath{IllumPara, Frame}{1});
    else
        %FP = Var.Analysis.FilePath{IllumPara,Frame}
        Img = imread(Var.Analysis.FilePath{IllumPara,Frame});
        [Var.Analysis.CurrentFolder, Var.Analysis.CurrentFileName]  = fileparts(Var.Analysis.FilePath{IllumPara, Frame});
    end
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(Img);
   % CurrentFN = Var.Analysis.CurrentFileName
    try
        Var.Analysis.TimeStamp = GetTimeStamp(Var);
    catch
        FileProps = dir(Var.Analysis.FilePath{IllumPara, Frame});
        TimeVec =   datevec(FileProps(1).datenum);
        Var.Analysis.TimeStamp = TimeVec(4)*3600 +TimeVec(5)*60+TimeVec(6) ;
    end
    
       % if Var.Analysis.CurrentFrame == 1 && CallNum == 1
%         fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
       try
            Var.Analysis.TimeZero = GetTimeStamp(Var, Var.Experiment.TimeZero);
        catch
            FileProps = dir(Var.Analysis.FilePath{IllumPara, Var.Experiment.TimeZero});
            TimeVecZero =   datevec(FileProps(1).datenum);
            Var.Analysis.TimeZero = TimeVecZero(4)*3600 +TimeVecZero(5)*60+TimeVecZero(6) ;
        end
   % end
    %% load from MicroManager from MultiDim from script data
elseif strcmp(Var.Experiment.Software, 'MicroManager') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimScript')
    Img = imread(Var.Analysis.FilePath{IllumPara,Frame});
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(Img);
    [Var.Analysis.CurrentFolder, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{IllumPara, Frame});
    % try
    Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = GetTimeStamp(Var);
    %     catch
    %       FileProps = dir(Var.Analysis.FilePath{IllumPara, Frame});
    %       TimeVec =   datevec(FileProps(1).datenum);
    %       Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = TimeVec(4)*3600 +TimeVec(5)*60+TimeVec(6) ;
    %     end
    
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
        % try
        Var.Analysis.TimeZero = GetTimeStamp(Var, Var.Experiment.TimeZero-Var.Analysis.FirstTPoint+1);
        %         catch
        %         FileProps = dir(Var.Analysis.FilePath{IllumPara, Var.Experiment.TimeZero});
        %         TimeVecZero =   datevec(FileProps(1).datenum);
        %         Var.Analysis.TimeZero = TimeVecZero(4)*3600 +TimeVecZero(5)*60+TimeVecZero(6) ;
        %         end
    end
    
    %% load from OPENLAB converted to Tiff by ImageJ
elseif strcmp(Var.Experiment.Software, 'OpenLab') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    FrameInStack = (Frame-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
    [TIFFImg] = tiffread_serge(Var.Analysis.FilePath{1}, FrameInStack);
    [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{1});
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
    end
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(TIFFImg.data);
    
    %% Load from TIFF stack
elseif strcmp(Var.Experiment.Software, 'TIF') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    %calculate the current Frame in the stack based on position
    %Illumination and frame number
    FrameInStack = (Frame-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
    %Read Tiff file
    [TIFFImg] = imread(Var.Analysis.FilePath{1}, FrameInStack);
    [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.Analysis.FilePath{1});
    if Frame == 1 && CallNum == 1
        fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
    end
    Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = double(TIFFImg);
    % Get time stamp from file info
    Var.Analysis.TimeStamp = GetTimeStamp(Var, FrameInStack);
    
    %Get Frame zero Time stamp
    FrameZeroInStack = (Var.Experiment.TimeZero-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
    Var.Analysis.TimeZero = GetTimeStamp(Var, FrameZeroInStack);
end

Var.Analysis.ImgSize = size(Var.Img.(Var.Analysis.LoadImgOut{CallNum}));


%     %Var.Analysis.TimeStamp(Var.Analysis.CurrentFrame) = TIFFImg.MM_stack(4);
% %% load from OpenLab software (Tiff files with numbers)
% elseif strcmp(Var.Microscope.Software, 'OL_num')
%     PosTime = (Frame-1)*(Var.ImgParameters.NumPos) + Pos;
%     [TIFFImg, img_read] = tiffread_serge(Var.ImgParameters.Filename{PosTime,IllumPara}, 1);
%     [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.ImgParameters.Filename{Pos,IllumPara});
%     if Frame == 1 && CallNum == 1
%         fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, '\n\n'])
%     end
% %% load from OpenLab software (Converted Liff)
% elseif strcmp(Var.Microscope.Software, 'OL')
%     FrameInStack = (Frame-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
%     [TIFFImg, img_read] = tiffread_serge(Var.ImgParameters.Filename{Pos,IllumPara}, FrameInStack);
%     [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.ImgParameters.Filename{Pos,IllumPara});
% %% load from tif stack
% elseif strcmp(Var.Microscope.Software, 'TIF') && strcmp(Var.Microscope.AcqType, 'SingleStack')
%     FrameInStack = (Frame-1)*(NumPos*NumIllum) + (Pos-1)*NumIllum + IllumPara;
%     [TIFFImg.data] = imread(Var.ImgParameters.Filename{Pos,IllumPara}, FrameInStack);
%     [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.ImgParameters.Filename{Pos,IllumPara});
%     if Frame == 1 && CallNum == 1
%         fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, ' at Position ' num2str(Pos) '\n\n'])
%     end
%
%     if length(size(TIFFImg.data)) >2
%         TIFFImg.data = TIFFImg.data(:,:,1);
%     end
% %% Load from Slidebook
% elseif strcmp(Var.Microscope.Software, 'SB')
%     %Check if data were acquired as a single file or multiple files for
%     %each time point
%     if strcmp(Var.ImgParameters.TimeLapseMode, 'Stack');
%         [TIFFImg, img_read] = tiffread_serge(Var.ImgParameters.Filename{Pos,IllumPara}, Frame);
%         [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.ImgParameters.Filename{Pos,1});
%         if Frame == 1 && CallNum == 1
%             fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, '\n\n'])
%         end
%     else
%         %Calculate filenumber for each file
%         FileNum = (Pos-1)*(Var.ImgParameters.NbFrame) + Frame;
%         %LoadedFile = Var.ImgParameters.Filename{FileNum,IllumSeg};
%         [TIFFImg, img_read] = tiffread_serge(Var.ImgParameters.Filename{FileNum,IllumSeg}, 1);
%         [pathstr, Var.Analysis.CurrentFileName] = fileparts(Var.ImgParameters.Filename{FileNum,1});
%         if CallNum == 1
%             fprintf(['\nAnalysing ', Var.Analysis.CurrentFileName, '\n\n'])
%         end
%     end


%% record image in variable

% if iscell(Var.Analysis.LoadImgOut{CallNum})
%     for Z = 1:length(Var.Analysis.LoadImageOut{CallNum})
%         Var.Img.(Var.Analysis.LoadImgOut{CallNum}{Z}) = double(TIFFImg.data{Z});
%     end
%     if strcmp(Var.Figure.Display, 'on')
%         FigNum = find(strcmp(Var.Figure.List, 'LoadImg'));
%         figure(FigNum(CallNum))
%         imagesc(Var.Img.(Var.Analysis.LoadImgOut{CallNum}{Z})); title (['loaded ', Var.Analysis.LoadIllum{CallNum}, ' image'])
%     end
% else


%% Display image
if strcmp(Var.Figure.Display, 'on')
    FigNum = find(strcmp(Var.Figure.List, 'LoadImg'));
    figure(FigNum(CallNum))
    %Set Colormap for figure
    if strcmp(Var.Experiment.Color{IllumPara}, 'B')
        CMap = [zeros(1,256); zeros(1,256); linspace(0,1,256)]';
        colormap(CMap)
    elseif strcmp(Var.Experiment.Color{IllumPara}, 'G')
        CMap = [zeros(1,256); linspace(0,1,256); zeros(1,256)]';
        colormap(CMap)
    elseif strcmp(Var.Experiment.Color{IllumPara}, 'R')
        CMap = [ linspace(0,1,256);zeros(1,256); zeros(1,256)]';
        colormap(CMap)
    else
        colormap(gray)
    end
    imagesc(Var.Img.(Var.Analysis.LoadImgOut{CallNum})); title (['loaded ', Var.Analysis.LoadIllum{CallNum}, ' image'])
    set(gca,'ytick',[], 'xtick', []);
end

else
    %If no new image loaded set NewImage to 0
    Var.Analysis.NewImage(IllumPara) = 0;
    %If no new image is loaded set it to an image filled with zeros
    if Var.Analysis.SkipFrame(1,IllumPara)<0
        Var.Img.(Var.Analysis.LoadImgOut{CallNum}) = zeros(Var.Analysis.ImgSize);
    end
end

Var.Analysis.SkipFrame(2,IllumPara) =  Var.Analysis.SkipFrame(2,IllumPara) + 1;
if Var.Analysis.SkipFrame(2,IllumPara) > abs(Var.Analysis.SkipFrame(1,IllumPara));
    Var.Analysis.SkipFrame(2,IllumPara) = 1;
end

%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;

%% SUB Function %%%%%%%%%%%%%%%

%% Get time stamp from MicroManager metadata

function TimeStamp = GetTimeStamp(Var, FrameNum)

if strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
    
    %get folder and file names
    %[Folder, FileName] = fileparts(ImgFilePath);
    Folder = Var.Analysis.CurrentFolder;
    if nargin == 1
        FileName = Var.Analysis.CurrentFileName;
    else
        if iscell(Var.Analysis.FilePath{1, FrameNum})
            [pathstr, FileName ] = fileparts(Var.Analysis.FilePath{1, FrameNum}{1});
        else
            [pathstr, FileName ] = fileparts(Var.Analysis.FilePath{1, FrameNum});
        end
    end
    
    
    %Open and read content of metadat file
    fid = fopen(fullfile(Folder,'metadata.txt'));
    C = textscan(fid, '%s','delimiter', '},');
    fclose(fid);
  %  assignin('base','Meta',C);
    %Find location of file name in metadata
    SearchStr = ['"FileName": "', FileName, '.tif"'];
    FindFileName = find(strcmp(C{1}, SearchStr));
    
    %find neighboring FrameKey
    %Look Down
    KeyFind = [];
    LookDown = FindFileName;
    while LookDown+1<length(C{1}) && isempty(KeyFind)
        LookDown = LookDown+1;
        KeyFind = cell2mat(strfind(C{1}(LookDown), 'FrameKey'));
    end
    %Look Up
    KeyFind = [];
    LookUp = FindFileName;
    while isempty(KeyFind)
        LookUp = LookUp-1;
        KeyFind = cell2mat(strfind(C{1}(LookUp), 'FrameKey'));
        
    end
    
    TimeInfo = [];
    while isempty(TimeInfo) && LookUp < LookDown
        LookUp = LookUp+1;
        TimeInfo = cell2mat(strfind(C{1}(LookUp), 'ElapsedTime'));
    end
    
    
    %Get time stamp
    NumLoc = isstrprop(C{1}{LookUp}, 'digit');
    TimeStamp = str2double(C{1}{LookUp}(NumLoc==1));
    
elseif strcmp(Var.Experiment.AcquisitionType, 'MultiDimScript')
    if nargin == 1
        FrameNum = Var.Analysis.CurrentFrame;
    end
    Folder = Var.Analysis.CurrentFolder;
    BaseDir = cd;
    %Get summary txt file
    cd(Folder)
    SumFile = dir('*_Summary.txt');
    if isempty(SumFile)
        cd ..
        SumFile = dir('*_Summary.txt');
    end
    
    fid = fopen(SumFile(1).name);
    C = textscan(fid, '%s','delimiter', '},');
    fclose(fid);
    
    cd(BaseDir);
    
    %Find location of file name in metadata
    SearchStr = ['Time ',num2str(FrameNum+Var.Analysis.FirstTPoint-1), ':'];
    FindTime = find(~cellfun('isempty', strfind(C{1}, SearchStr)));
    %Get String
    TimeStr = C{1}{FindTime};
    Col = strfind(TimeStr, ':');
    %Get time value
    TimeStamp = str2double(TimeStr(Col+1:end));
elseif strcmp(Var.Experiment.Software, 'TIF') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    %Get TIFF File info
    Info = imfinfo(Var.Analysis.FilePath{1});
    if isfield(Info, 'ImageDescription')
        ImageMetaData = Info(FrameNum).ImageDescription;
        if strcmp(ImageMetaData(2:17), 'Improvision Data')
            TSstr = strfind(ImageMetaData, 'TimeStamp');
            TimeStamp = str2double(ImageMetaData(TSstr(1)+22:TSstr(2)-2));
        end
    end
    
end



%% 

%%%%%%%%%%%%%%%%%%%%%

%% LOAD Z STACK FCT

function MaxImg = LoadZstack(FileList)
%Load img from file names given in cell and makes maximum intensity
%projection
for Z = 1:length(FileList)
    Img = imread(FileList{Z});
    if Z == 1
        MaxImg = Img;
    else
        MaxImg = max(MaxImg,Img);
    end
%                         figure(100)
%                         imagesc(MaxImg); title(num2str(Z));
%                         colormap(gray)
%                         drawnow
end %end Z scan
