function [Var] = ReadImgParameters(Var)

%% Go to Img folder and find files
%Save current dir as base directory
BaseDir = cd;
%Get Folder name from Experiment structure
ImgFolder = Var.Experiment.Folder;
ImgFolderTxt = regexprep(ImgFolder, '\', '\\\\');
cd(ImgFolder)

%get Image name
ImgName = Var.Experiment.Filename;

%Get position and Illum
NumPos = Var.Analysis.NumPos;
NumPosAnalysis = Var.Analysis.NumPosAnalyzed;
Var.Analysis.NumIllum = length(Var.Experiment.Illumination);
NumIllum = Var.Analysis.NumIllum;

if isfield(Var.Experiment, 'FinalTime')
    MaxTimePoint = Var.Experiment.FinalTime;
else
    MaxTimePoint = Inf;
end

%set extension based  software and acquisition type

if  strcmp(Var.Experiment.Software, 'MetaMorph')
    if strcmp(Var.Experiment.AcquisitionType, 'IterSave') || strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
        Extension = '*stk';
        SearchString = [ImgName, Extension];
    elseif strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
        Extension = '*TIF';
        SearchString = [ImgName, Extension];
    end
elseif strcmp(Var.Experiment.Software, 'MicroManager')
    if strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional') || strcmp(Var.Experiment.AcquisitionType, 'MultiDimScript')
        SearchString = [ImgName, '*tif'];
        TifFiles = dir(SearchString);
        
        if ~isempty(TifFiles)
            PositionFolderName = [];
        else
            PositionDir = dir(ImgFolder);
            
            %assignin('base', 'Pdir', PositionDir)
            
            IterPos = 0;
            for P = 1:length(PositionDir)
                
                PosStr = strfind(PositionDir(P).name, 'Pos');
                
                if PositionDir(P).isdir && ~isempty(PosStr) && PosStr == 1
                    
                    PosNumStr = [];
                    PosStr = PosStr+3;
                    % PD = PositionDir(P).name
                    %Dig = isstrprop(PositionDir(P).name(PosStr), 'digit')
                    while PosStr<= length(PositionDir(P).name) && isstrprop(PositionDir(P).name(PosStr), 'digit');
                        PosNumStr = [PosNumStr,PositionDir(P).name(PosStr)];
                        PosStr = PosStr+1;
                    end
                    %PosFolderNum = str2double(PositionDir(P).name(NumLoc==1));
                    PosFolderNum = str2double(PosNumStr);
                    %check for the first folder if the Pos numbering start from 0 or 1
                    if IterPos == 0
                        if PosFolderNum == 0
                            AddOne = 1;
                            IterPos = 1;
                        else
                            AddOne = 0;
                            IterPos = 1;
                        end
                    end
                    
                    if Var.Analysis.CurrentPos == PosFolderNum+AddOne  %+1 if pos number start with 0
                        SearchString = [PositionDir(P).name, filesep, ImgName, '*tif'];
                        PositionFolderName = PositionDir(P).name;
                        IterPos = IterPos+1;
                        %                                     CP = Var.Analysis.CurrentPos
                        %                 PosF = PosFolderNum+AddOne
                    end
                end
            end
        end
        %Check if Z-stack is acquired and for which illum
        Var.Analysis.Zstack = ones(1,Var.Analysis.NumIllum);
        if isfield(Var.Experiment, 'Zstart')
            LengthZstart = length(Var.Experiment.Zstart);
            for Z = 1:LengthZstart
                if ~isempty(Var.Experiment.Zstart{Z})
                    %Calculate number of planes
                    Var.Analysis.Zstack(Z) = 1+floor((Var.Experiment.Zend{Z}-Var.Experiment.Zstart{Z})/Var.Experiment.Zstep{Z});
                end
            end
        end
        
        %Check if Frame skipped and for which illum
        Var.Analysis.SkipFrame = ones(2,Var.Analysis.NumIllum);
        if isfield(Var.Experiment, 'Skip')
            LengthSkip = length(Var.Experiment.Skip);
            for S = 1:LengthSkip
                if ~isempty(Var.Experiment.Skip{S})
                    Var.Analysis.SkipFrame(1,S) = Var.Experiment.Skip{S}+1;
                    if isfield(Var.Experiment, 'SkipSA') && strcmp(Var.Experiment.SkipSA{S}, 'Analyze')
                        Var.Analysis.SkipFrame(1,S) = -1*Var.Analysis.SkipFrame(1,S);
                    end
                        
                end
            end
        end
        
                
        
        %Check file list for lowest timepoint saved
        if strcmp(Var.Experiment.AcquisitionType, 'MultiDimScript')
            ImgFiles = dir(SearchString);
            for f = 1:length(ImgFiles)
                FileName = ImgFiles(f).name;
                %dissect filename to extract illumination and Time point
                US_T = strfind(FileName, '_T');
                TimePoint = str2double(FileName(US_T(end)+2:US_T(end)+4));
                AllT(f) = double(TimePoint);
            end
            %Set lowest time point to FirstTpoint
            Var.Analysis.FirstTPoint = min(AllT);
        end
    end
    
elseif  strcmp(Var.Experiment.Software, 'NikonElement') && (strcmp(Var.Experiment.AcquisitionType, 'TiffSerie') || strcmp(Var.Experiment.AcquisitionType, 'TiffExport'))
    Extension = '*tif';
    SearchString = [ImgName, Extension];
elseif  strcmp(Var.Experiment.Software, 'OpenLab') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    Extension = '*tif';
    SearchString = [ImgName, Extension];
    
elseif  strcmp(Var.Experiment.Software, 'TIF') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    Extension = '*tif';
    SearchString = [ImgName, Extension];
     Var.Analysis.SkipFrame = ones(2,Var.Analysis.NumIllum);
end

% Find files with matching initial string in ImgFolder
ImgFiles = dir(SearchString);

%if no file is found
if isempty(ImgFiles)
    error(['No files found looking for ' SearchString ' in ' ImgFolderTxt])
end

cd(BaseDir)

assignin('base', 'ImgFiles', ImgFiles)

%% Distrubute files for Image analysis
NbFiles = size(ImgFiles,1);



% get files depending on different software and acquisition types
if  strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    % For single stack one should have a single image file
    if length(ImgFiles) > 1
        error(['Multiple files found with ' SearchString ' in ' ImgFolderTxt])
    end
    %Generate full filepath
    Var.Analysis.FilePath{1} = fullfile(ImgFolder, ImgFiles.name);
    %Generate path for saving analysis
    Dot = strfind(ImgFiles.name, '.');
    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles.name(1:Dot-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
    
    %Try to read images
    [stack, img_read] = tiffread_serge(Var.Analysis.FilePath{1}, 0);
    Var.Analysis.NumFrame = img_read/(NumPos*NumIllum);
    
    
elseif  strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'IterSave')
    %For IterSave multiple file will have the same root Use first one for as
    %basename but remove number
    US  = strfind(ImgFiles(1).name,'_');
    for i = 1:NbFiles
        %Generate full filepath
        Var.Analysis.FilePath{i} = fullfile(ImgFolder, [ImgFiles(1).name(1:US(end)), num2str(i), '.stk']);
    end
    %Generate path for saving analysis
    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles(1).name(1:US(end)-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
    %Get Number of frames for number of files found
    Var.Analysis.NumFrame = NbFiles;
    
elseif  strcmp(Var.Experiment.Software, 'MetaMorph') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
    AllPos = [];
    NbFilesAnalysis = 0;
    for f = 1:NbFiles
        %Get Filename
        FileName = ImgFiles(f).name;
        %dissect filename to extract Position, illumination, Time point
        US = strfind(FileName, '_');
        Dot = strfind(FileName, '.');
        TimePoint = str2double(FileName(US(end)+2:Dot(end)-1));
        Position = str2double(FileName(US(end-1)+2:US(end)-1));
        Illum = FileName(US(end-2)+1:US(end-1)-1);
        
        %Check if file position matches current position
        if Var.Analysis.CurrentPos == Position
            for I = 1:NumIllum
                if strfind(Illum, Var.Experiment.Filter{I})
                    T = double(TimePoint);
                    Var.Analysis.FilePath{I,T} = fullfile(ImgFolder, ImgFiles(f).name);
                    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles(f).name(1:US(end-2)-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
                    NbFilesAnalysis = NbFilesAnalysis +1;
                end
            end
        end
        %if not all positions are present for analysis
        AllPos = [AllPos, Position];
        
        
        
    end
    %Remove duplicated position numbers from AllPos
    Var.Analysis.AllPos = unique(AllPos);
    
    %Get true number of positions in folder
    NumPos = length(Var.Analysis.AllPos);
    
    %get number of frames from MD data
    Var.Analysis.NumFrame = NbFilesAnalysis/(NumPos*NumIllum);
    %% MICROMAMAGER  MULTI Dim acq
elseif  strcmp(Var.Experiment.Software, 'MicroManager') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimensional')
    NbFilesAnalysis = 0;
    for f = 1:NbFiles
        %Get Filename
        FileName = ImgFiles(f).name;
        %dissect filename to extract illumination and Time point
        US = strfind(FileName, '_');
        if f == 1
            Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles(f).name(1:US(end-2)-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
        end
        
        %Get time point and illum name
        TimePoint = str2double(FileName(US(1)+1:US(2)-1))+1;
        Illum = FileName(US(2)+3:US(end));
        
        %Check which illum is in filename
        for I = 1:NumIllum
            if strfind(Illum, Var.Experiment.Filter{I})
                T = double(TimePoint);
                %Test if time point is below max time point
                if T <= MaxTimePoint
                %Check if single Z plane
                if Var.Analysis.Zstack(I) == 1
                    Var.Analysis.FilePath{I,T} = fullfile(ImgFolder, PositionFolderName, ImgFiles(f).name);
                    NbFilesAnalysis = NbFilesAnalysis +1;
                else %For multiple Z plane put all files in a cell
                    %Get slice number from filename
                    Zslice = str2double(FileName(US(3)+1:end-4))+1;
                    Var.Analysis.FilePath{I,T}{Zslice} = fullfile(ImgFolder, PositionFolderName, ImgFiles(f).name);
                    NbFilesAnalysis = NbFilesAnalysis +1;
                end
                end
            end
            
        end
    end
    
    %get number of frames from MD data
    %Var.Analysis.NumFrame = NbFilesAnalysis/NumIllum;
    Var.Analysis.NumFrame = NbFilesAnalysis/sum(Var.Analysis.Zstack);

    
    
    
    %verify with size FilePath Cell
    NumFilePath = size( Var.Analysis.FilePath,2);
    
    if Var.Analysis.NumFrame ~= NumFilePath
        WarnStr = ['Mismatch between Number of frame calculated: ', num2str(Var.Analysis.NumFrame), ...
            ' and Number of files in FilePath matrix: ',num2str(NumFilePath)];
        warning(WarnStr)
        
        Var.Analysis.NumFrame = NumFilePath;
    end
    
%%  MICROMANAGER Multidim script
elseif  strcmp(Var.Experiment.Software, 'MicroManager') && strcmp(Var.Experiment.AcquisitionType, 'MultiDimScript')
    NbFilesAnalysis = 0;
    for f = 1:NbFiles
        %Get Filename
        FileName = ImgFiles(f).name;
        %dissect filename to extract illumination and Time point
        US_T = strfind(FileName, '_T');
        TimePoint = str2double(FileName(US_T(end)+2:US_T(end)+4));
        IllumStart = length(ImgName)+2;
        Illum = FileName(IllumStart:US_T(end)-1);
        
        %Check if file position matches current position
        
        for I = 1:NumIllum
            if strfind(Illum, Var.Experiment.Filter{I})
                T = double(TimePoint) -(Var.Analysis.FirstTPoint)+1;
                Var.Analysis.FilePath{I,T} = fullfile(ImgFolder, PositionFolderName, ImgFiles(f).name);
                Var.Analysis.OutPath = [fullfile(ImgFolder, ImgName), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
                NbFilesAnalysis = NbFilesAnalysis +1;
            end
        end
        
    end
    
    %get number of frames from MD data
    Var.Analysis.NumFrame = NbFilesAnalysis/NumIllum;
    
elseif  strcmp(Var.Experiment.Software, 'NikonElement') && strcmp(Var.Experiment.AcquisitionType, 'TiffSerie')
    
    for f = 1:NbFiles
        %Get Filename
        FileName = ImgFiles(f).name;
        %dissect filename to extract Position, illumination, Time point
        T = strfind(FileName, 'timet');
        Dot = strfind(FileName, '.');
        TimePoint = str2double(FileName(T(end)+6:Dot(end)-1));
        
        
        T = double(TimePoint);
        Var.Analysis.FilePath{1,T} = fullfile(ImgFolder, ImgFiles(f).name);
        Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles(f).name(1:T(end)-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
    end
    
    %get number of frames from data
    Var.Analysis.NumFrame = NbFiles;
    
    
elseif  strcmp(Var.Experiment.Software, 'NikonElement') && strcmp(Var.Experiment.AcquisitionType, 'TiffExport')
    AllPos = [];
    for f = 1:NbFiles
        %Get Filename
        FileName = ImgFiles(f).name;
        %dissect filename to extract Position, illumination, Time point
        US = strfind(FileName, '_p');
        t0 = strfind(FileName, 't0');
        Dot = strfind(FileName, '.');
        
        Position = str2double(FileName(US(end)+2:US(end)+3));
        
        TimePoint = str2double(FileName(t0+1:t0+5));
        
        Illum = FileName(Dot(end)-3:Dot(end)-1);
        
        %Check if file position matches current position
        if Var.Analysis.CurrentPos == Position
            for I = 1:NumIllum
                if strfind(Illum, Var.Experiment.Filter{I})
                    T = double(TimePoint);
                    Var.Analysis.FilePath{I,T} = fullfile(ImgFolder, ImgFiles(f).name);
                    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles(f).name(1:US(end)-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
                end
            end
        end
        %if not all positions are present for analysis
        AllPos = [AllPos, Position];
        
        
        
    end
    %Remove duplicated position numbers from AllPos
    Var.Analysis.AllPos = unique(AllPos);
    
    %Get true number of positions in folder
    NumPos = length(Var.Analysis.AllPos);
    
    %get number of frames from MD data
    Var.Analysis.NumFrame = NbFiles/(NumPos*NumIllum);
elseif  strcmp(Var.Experiment.Software, 'TIF') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    % For single stack one should have a single image file
    if length(ImgFiles) > 1
        error(['Multiple files found with ' SearchString ' in ' ImgFolderTxt])
    end
    NbFilesAnalysis = 1;
    %Generate full filepath
    Var.Analysis.FilePath{1} = fullfile(ImgFolder, ImgFiles.name);
    %Generate path for saving analysis
    Dot = strfind(ImgFiles.name, '.');
    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles.name(1:Dot-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
    
    %Try to read images
    ImgInfo = imfinfo(Var.Analysis.FilePath{1});
    Var.Analysis.NumFrame = numel(ImgInfo)/(NumPos*NumIllum);
    
elseif  strcmp(Var.Experiment.Software, 'OpenLab') && strcmp(Var.Experiment.AcquisitionType, 'SingleStack')
    % For single stack one should have a single image file
    if length(ImgFiles) > 1
        error(['Multiple files found with ' SearchString ' in ' ImgFolderTxt])
    end
    NbFilesAnalysis = 1;
    %Generate full filepath
    Var.Analysis.FilePath{1} = fullfile(ImgFolder, ImgFiles.name);
    %Generate path for saving analysis
    Dot = strfind(ImgFiles.name, '.');
    Var.Analysis.OutPath = [fullfile(ImgFolder, ImgFiles.name(1:Dot-1)), '_Pos', num2str(Var.Analysis.CurrentPos,'%02d')];
    
    %Try to read images
    [stack, img_read] = tiffread_serge(Var.Analysis.FilePath{1}, 0);
    Var.Analysis.NumFrame = img_read/(NumPos*NumIllum);
    
end



%Output number files found and number illumination and number of frame
if Var.Analysis.CurrentPos == NumPos
    
    
    if strcmp(Var.Experiment.Software, 'MicroManager')
        
        %ResultStr  = [num2str(NbFilesAnalysis), ' Files named ', SearchString ' found in: ', ImgFolderTxt, ' reprensenting:\n' num2str(NumPos) ' Positions ', num2str(NumIllum), ' Exposure types and ' , num2str(Var.Analysis.NumFrame), ' Frames\n\n'];
        
        ResultStr = ['Found ',  num2str(NbFilesAnalysis), ' Files in folder ', ImgFolderTxt, ' searching for ', SearchString ,'\n',...
            'Representing: 1 of ' ,num2str(NumPos), ' Stage Positions for ' , num2str(Var.Analysis.NumFrame), ' Time Frames and ', num2str(NumIllum), ' Exposure types\n'];
        if isfield(Var.Analysis, 'Zstack') && max(Var.Analysis.Zstack) >1
            for I = 1:NumIllum
                ResultStr = [ResultStr, Var.Experiment.Filter{I}, ': ', num2str(Var.Analysis.Zstack(I)), ' Z-slice\n'];
            end
        end
        
        fprintf(ResultStr)
    else
        ResultStr  = [num2str(NbFilesAnalysis), ' Files named ', SearchString ' found in: ', ImgFolderTxt, ' reprensenting:\n' num2str(NumPos) ' Positions ', num2str(NumIllum), ' Exposure types and ' , num2str(Var.Analysis.NumFrame), ' Frames\n\n'];
        fprintf(ResultStr)
    end
end

