function VarCell = PrepareAnalysis(SearchField, System, PosPerWell)
%% Function to get experiments parameters from the database and prepare the files for the analysis
%SearchField:   Either a vector containing the experiment Numbers
%               Or a Cell array where the first item is the Search criteria and the second the Search string
%SearchStr: String that contains characters to search for
%  or cell that contains multiple strings to search for
%System :   mac Run on own machine
%           pc run on own machine
%           Linux: Run on signaling server
%           HPC: Run on argos
%           Wally: run on Wally
%PosPerWell: Sets the number of position that are analyzed for each well.

%Example:VarCell = PrepareAnalysis( {'Filename', '0326_Hog1_WT'}, 'mac')
%		 VarCell = PrepareAnalysis( [22:25], 'Linux')

%% Initiate routine by setting directory, user specific variables and getting required inputs

UserName = 'spelet1';

% Set of File paths: Might change depending on user.
WallyPath = '/scratch/wally/FAC/FBM/DMF/spelet1';
HPCPath = '/groups/dmf/gr_Pelet';   
LinuxPath = ['/mnt/',UserName, '/'];
CellXcore = 'CellX/core';
try
    %Add path to folder containing all routines
    BaseDir = cd;
    addpath(fullfile(BaseDir, 'CommonProcess'))
    %Load database connection settings from install process
    Var = load('DatabaseConnextionParameter.mat');
    %DBpath = Var.Database.DriverPath
    javaaddpath(Var.Database.DriverPath)
    load('ImageStoragePath.mat');
catch
    error('Could not find DatabaseConnextionParameter.mat or ImageStoragePath.mat. Verify that these files are in the CommonProcess folder')
end

%% Dialog for getting the parameters
%If search argument for database not provided display dialog box
if nargin == 0
    DialogPrompt = {'Enter Search Field (Folder, Filename):','Enter Search String:'};
    DialogTitle = 'Input for Database search';
    DialogDef = {'Filename','{0326_Hog1_WT}'};
    DialogAnswer = inputdlg(DialogPrompt, DialogTitle,1 ,DialogDef);
    
    Var.Database.SearchField = DialogAnswer{1};
    Var.Database.SearchStr = DialogAnswer{2};
else
    %Check if only a vector is provided => Search criteria is ExpNum
    if isnumeric(SearchField)
        Var.Database.SearchField = 'ExpNum';
        Var.Database.SearchStr = SearchField;
    else
        Var.Database.SearchField = SearchField{1};
        Var.Database.SearchStr = SearchField{2};
    end
end

%Select system to run the analysis
if nargin < 2
    System = questdlg('Select system to run the analysis','System choice','mac','Linux','HPC', 'Wally', 'mac') ;
end

if nargin < 3
    PosPerWell = Inf;
end

try
    
    %% Get Experiments Data from database connection
    Var.Database.Table = 'Experiment';
    Var = DatabaseConnections(Var);
    assignin('base', 'Var_exp', Var)
    
    %% Iteration variable to build the cell
    iter = 0;
    for R = 1:Var.Database.NumFoundRecord
        %Distribute database fields into VarRec
        VarRec = Var;
        VarRec = GetDBFields(VarRec, R);
        %% Read Analysis info from database
        VarRec.Database.Table = 'Analysis';
        VarRec.Database.SearchField = 'AnalysisNum';
        VarRec.Database.SearchStr = VarRec.Experiment.LinkedAnalysisNum;
        %Get analysis flow parameters from database
        VarRec = DatabaseConnections(VarRec);
        
        if VarRec.Database.NumFoundRecord ~= 1
            error(['Found ' , VarRec.Database.NumFoundRecord, ' for Analysis Number ' , VarRec.Experiment.LinkedAnalysisNum,' Please check the Database'])
        end
        
        %Distribute database fields into VarRec
        VarRec = GetDBFields(VarRec, 1);
        
        VarRec.Analysis.System = System;
        VarRec.Analysis.RunningUser = UserName;
        
        %Get server folder path for window and Mac platform
        %from the ImageStoragePath file
        %Use local if you run the process on the local machine
        
        VarRec.Analysis.MacServerBase =  ImageStoragePath.MAC;  %'local';  %
        VarRec.Analysis.PCServerBase = ImageStoragePath.PC;        % 'local';  %
        
        
        %For system running in parallel no movies can be saved jpg of each
        %frame will be saved instead
        VarRec.Analysis.SaveAvi = 'No';
        %     if ~isempty(strfind(VarRec.Analysis.System, 'para'))
        %         VarRec.Analysis.SaveAvi = 'No';
        %     else
        %         VarRec.Analysis.SaveAvi = 'Yes';
        %     end
        
        
        %Find Pathformat of DB Folder entry
        if strcmp(VarRec.Analysis.MacServerBase, 'local')
            %Keep current path if no server path is specified
            VarRec.Analysis.MacPath = VarRec.Experiment.Folder;
            VarRec.Analysis.PCPath = VarRec.Experiment.Folder;
        else
            FindMacSep = strfind(VarRec.Experiment.Folder, '/');
            if ~isempty(FindMacSep)
                VarRec.Analysis.MacPath = VarRec.Experiment.Folder;
              %  VarRec.Analysis.PCPath = Mac2PC(VarRec.Experiment.Folder, VarRec.Analysis.PCServerBase);
            else
                VarRec.Analysis.PCPath = VarRec.Experiment.Folder;
                VarRec.Analysis.MacPath = PC2MAC(VarRec.Experiment.Folder, VarRec.Analysis.MacServerBase);
            end
        end
        
        %assignin('base', 'Var_2', VarRec)
        
        %% Read segmentation settings from database
        VarRec.Database.Table = 'SegmentationParameters';
        VarRec.Database.SearchField = 'FullName';
        if isfield(VarRec.Analysis, 'SegYFPara')
            % ShStr = [VarRec.Analysis.SegYFPara{1},'_',num2str(VarRec.Experiment.Objective),'x_bin', num2str(VarRec.Experiment.Bin)]
            VarRec.Database.SearchStr = [VarRec.Analysis.SegYFPara{1},'_',num2str(VarRec.Experiment.Objective),...
                'x_bin', num2str(VarRec.Experiment.Bin)];
            VarRec = DatabaseConnections(VarRec);
            %Distribute database fields into VarRec
            VarRec = GetDBFields(VarRec, 1);
        end
        if isfield(VarRec.Analysis, 'SegAroundPara')
            VarRec.Database.SearchStr = [VarRec.Analysis.SegAroundPara{1},'_',num2str(VarRec.Experiment.Objective),...
                'x_bin', num2str(VarRec.Experiment.Bin)];
            VarRec = DatabaseConnections(VarRec);
            %Distribute database fields into VarRec
            VarRec = GetDBFields(VarRec, 1);
        end
        if isfield(VarRec.Analysis, 'SegPhasePara')
            VarRec.Database.SearchStr = [VarRec.Analysis.SegPhasePara{1},'_',num2str(VarRec.Experiment.Objective),...
                'x_bin', num2str(VarRec.Experiment.Bin)];
            VarRec = DatabaseConnections(VarRec);
            %Distribute database fields into VarRec
            VarRec = GetDBFields(VarRec, 1);
        end
        
        if isfield(VarRec.Analysis, 'SegRNAPara')
            VarRec.Database.SearchStr = [VarRec.Analysis.SegRNAPara{1},'_',num2str(VarRec.Experiment.Objective),...
                'x_bin', num2str(VarRec.Experiment.Bin)];
            VarRec = DatabaseConnections(VarRec);
            %Distribute database fields into VarRec
            VarRec = GetDBFields(VarRec, 1);
        end
        
        if strcmp(VarRec.Experiment.Segmentation, 'CellXNucl')
            %Load CellX Object
            addpath(genpath(CellXcore));
            load('CellXObj.mat')
            VarRec.CellXObj = CellXObj;
        end
        
        %% Assign figure numbers depending on analysis type
        VarRec = FigureFlow(VarRec);
        
        %% Creat Cell array containing data for each parrallel analysis
        %Get Number of position from Experiment parameters
        VarRec.Analysis.NumPos = VarRec.Experiment.LastPos{length(VarRec.Experiment.LastPos)};
        VarRec.Analysis.NumPosAnalyzed = min(VarRec.Analysis.NumPos, length(VarRec.Experiment.Well)*PosPerWell);
        %Get Well number and condition from experiment data
        
        
        LineNum = 1;
        %Loop Through all positions
        for P = 1:VarRec.Analysis.NumPos
            %P = P
            %Set Current position number to P
            VarRec.Analysis.CurrentPos = P;
            %Check if P is within the position bounds of the current line
            if P >= VarRec.Experiment.FirstPos{LineNum} && P <= VarRec.Experiment.LastPos{LineNum}
                %Check if P is within the number of positions to be
                %analyzed in one well
                if P >= VarRec.Experiment.FirstPos{LineNum} && P <= VarRec.Experiment.FirstPos{LineNum}+PosPerWell-1
                    VarRec.Analysis.CurrentWell =  VarRec.Experiment.Well{LineNum};
                    VarRec.Analysis.CurrentCondition =  VarRec.Experiment.Condition{LineNum};
                else
                    %If it is larger: skip remaining part of the loop
                    fprintf(['Skip Pos ', num2str(P), '\n'])
                    continue
                end
            else
                %Add one to line number to move to next line in database
                LineNum = LineNum+1;
                if P >= VarRec.Experiment.FirstPos{LineNum} && P <= VarRec.Experiment.LastPos{LineNum}
                    if P >= VarRec.Experiment.FirstPos{LineNum} && P <= VarRec.Experiment.FirstPos{LineNum}+PosPerWell-1
                        VarRec.Analysis.CurrentWell =  VarRec.Experiment.Well{LineNum};
                        VarRec.Analysis.CurrentCondition =  VarRec.Experiment.Condition{LineNum};
                    else
                        fprintf(['Skip Pos ',  num2str(P), '\n'])
                        continue
                    end
                else
                    error('Position number not found in position range from experiment')
                end
            end
            
            
            
            %Read Image variable from files and define image path and output
            %path
            %Change the folder Path depending on platform where the analysis preparation is performed
            if ismac
                VarRec.Experiment.Folder = VarRec.Analysis.MacPath;
            elseif ispc
                VarRec.Experiment.Folder = VarRec.Analysis.PCPath;
            end
            %Populate VarCell with current Variable for the record
            iter = iter + 1;
            VarCell{iter} = VarRec;
            
            %Get Image parameter from image folder and add to Cell
            [VarCell{iter}] = ReadImgParameters(VarCell{iter});
            %assignin('base', 'Var_3', VarRec)
        end
        
    end
catch LastErr
    VarRec.error = LastErr;
    assignin('base','VarRec_Error',VarRec);
    %Save VarRec file
    Time = clock;
    TimeStr = [num2str(Time(1)),num2str(Time(2),'%02d'),num2str(Time(3),'%02d'),'_',num2str(Time(4),'%02d'),num2str(Time(5),'%02d')];
    save(fullfile(BaseDir,['ERROR_',TimeStr,'.mat']), 'VarRec')
    
    
end
assignin('base', 'VarCell_RIP', VarCell)

%% Flatness analysis
%VarCell = FlatnessMeasure(VarCell);

%% Save Cell array depending on analysis platform used

if strcmpi(System, 'Wally')
    %Get Number of positions
    NumPositions = length(VarCell);
    %Calculate the number of frames to be segmented
    TotalNumFrame = 0;
    for P = 1:NumPositions
        if isfield(VarCell{P}.Experiment, 'Skip')
            VarCell{P}.Analysis.NumSkip = max(VarCell{P}.Analysis.SkipFrame(1,:));
            
            TotalNumFrame = TotalNumFrame + ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
            VarCell{P}.Analysis.NumFrame2Analyze = ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
        else
            TotalNumFrame = TotalNumFrame + VarCell{P}.Analysis.NumFrame;
            VarCell{P}.Analysis.NumFrame2Analyze = VarCell{P}.Analysis.NumFrame;
            VarCell{P}.Analysis.NumSkip = 1;
        end
    end
    % assignin('base', 'VarCell_RIP', VarCell)
    %Loop through all frames
    for F = 1:TotalNumFrame
        
        % Calculate Position and Frame number
        CurrentFrame = F;
        CurrentCell = 1;
        while CurrentCell<NumPositions && CurrentFrame > VarCell{CurrentCell}.Analysis.NumFrame2Analyze
            CurrentFrame = CurrentFrame - VarCell{CurrentCell}.Analysis.NumFrame2Analyze;
            CurrentCell = CurrentCell+1;
        end
        
        %If Skipped frames:
        if VarCell{CurrentCell}.Analysis.NumSkip >1
            CurrentFrame = VarCell{CurrentCell}.Analysis.NumSkip*(CurrentFrame-1)+1;
        end
        
        %Transfer Position and Frame numbers to VarFrame
        VarFrame = VarCell{CurrentCell};
        VarFrame.Analysis.CurrentCell = CurrentCell;
        VarFrame.Analysis.CurrentFrame = CurrentFrame;
        %Check if Skipped frame for analysis
        
        if CurrentFrame > 1 && isfield(VarFrame.Experiment, 'Skip') && min(VarFrame.Analysis.SkipFrame(1,:))< 0
            
            SkipAnalysis = find(VarFrame.Analysis.SkipFrame(1,:)<0);
            for SA = 1:length(SkipAnalysis)
                VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) =  VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) + 1;
                if VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) > abs(VarFrame.Analysis.SkipFrame(1,SkipAnalysis(SA)));
                    VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) = 1;
                end
                VarCell{CurrentCell}.Analysis.SkipFrame(2,SkipAnalysis(SA)) = VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA));
            end
            
        end
        
        %Generate Var File name
        FrameName = [VarFrame.Experiment.Filename,'Pos', num2str(VarFrame.Analysis.CurrentPos,'%03.d_Frame'), num2str(VarFrame.Analysis.CurrentFrame, '%05d') ,'_VAR.mat'];
        
        %Modify path for Wally
        VarFrame = ChangPathForWally(VarFrame, WallyPath);
        %Generate HPC folder
        if ~exist(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC'), 'file')
            mkdir(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC'))
        end
        %Save VAR file to HPC folder
        save(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC',FrameName) , 'VarFrame');
        
        %Text file to write names of VarFile
        if F == 1
            %Create Time stamp
            Now = clock;
            TimeStamp = [num2str(Now(1)),num2str(Now(2),'%02.0f'),num2str(Now(3),'%02.0f'),'_',num2str(Now(4),'%02.0f'),num2str(Now(5),'%02.0f') ];
            
            %text File Name
            VarListFilepath = fullfile(VarCell{CurrentCell}.Experiment.Folder,[TimeStamp, '_VARList.txt']);
            %Open txt file
            fileID = fopen(VarListFilepath, 'a+','n','UTF-8');
            %Add Var name to txt
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
            
            %Write JobArray script file
            JobID = ['EN', num2str(VarFrame.Experiment.ExpNum), '_', VarFrame.Experiment.User];
            Wally_ScriptWrite(JobID, TotalNumFrame, fullfile(VarFrame.Analysis.MacPath,[TimeStamp, '_VARList.txt']), VarCell{CurrentCell}.Experiment.Folder)
            Wally_WrapperWrite(JobID, TimeStamp, VarCell)
            
        elseif F == TotalNumFrame
             %Add Var name to txt
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
            %Close text
            fclose(fileID);
            
        else
             %Add Var name to txt
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
        end
        
        
        
    end    
    
elseif strcmpi(System, 'HPC')
    %Get Number of positions
    NumPositions = length(VarCell);
    %Calculate the number of frames to be segmented
    TotalNumFrame = 0;
    for P = 1:NumPositions
        if isfield(VarCell{P}.Experiment, 'Skip')
            VarCell{P}.Analysis.NumSkip = max(VarCell{P}.Analysis.SkipFrame(1,:));
            
            TotalNumFrame = TotalNumFrame + ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
            VarCell{P}.Analysis.NumFrame2Analyze = ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
        else
            TotalNumFrame = TotalNumFrame + VarCell{P}.Analysis.NumFrame;
            VarCell{P}.Analysis.NumFrame2Analyze = VarCell{P}.Analysis.NumFrame;
            VarCell{P}.Analysis.NumSkip = 1;
        end
    end
    % assignin('base', 'VarCell_RIP', VarCell)
    %Loop through all frames
    for F = 1:TotalNumFrame
        
        % Calculate Position and Frame number
        CurrentFrame = F;
        CurrentCell = 1;
        while CurrentCell<NumPositions && CurrentFrame > VarCell{CurrentCell}.Analysis.NumFrame2Analyze
            CurrentFrame = CurrentFrame - VarCell{CurrentCell}.Analysis.NumFrame2Analyze;
            CurrentCell = CurrentCell+1;
        end
        
        %If Skipped frames:
        if VarCell{CurrentCell}.Analysis.NumSkip >1
            CurrentFrame = VarCell{CurrentCell}.Analysis.NumSkip*(CurrentFrame-1)+1;
        end
        
        %Transfer Position and Frame numbers to VarFrame
        VarFrame = VarCell{CurrentCell};
        VarFrame.Analysis.CurrentCell = CurrentCell;
        VarFrame.Analysis.CurrentFrame = CurrentFrame;
        %Check if Skipped frame for analysis
        
        if CurrentFrame > 1 && isfield(VarFrame.Experiment, 'Skip') && min(VarFrame.Analysis.SkipFrame(1,:))< 0
            
            SkipAnalysis = find(VarFrame.Analysis.SkipFrame(1,:)<0);
            for SA = 1:length(SkipAnalysis)
                VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) =  VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) + 1;
                if VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) > abs(VarFrame.Analysis.SkipFrame(1,SkipAnalysis(SA)));
                    VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA)) = 1;
                end
                %   assignin('base', ['VarF_',num2str(F)], VarFrame)
                VarCell{CurrentCell}.Analysis.SkipFrame(2,SkipAnalysis(SA)) = VarFrame.Analysis.SkipFrame(2,SkipAnalysis(SA));
            end
            
        end
        
        
        FrameName = [VarFrame.Experiment.Filename,'Pos', num2str(VarFrame.Analysis.CurrentPos,'%03.d_Frame'), num2str(VarFrame.Analysis.CurrentFrame, '%05d') ,'_VAR.mat'];
        
        VarFrame = ChangPathForCluster(VarFrame, HPCPath);
        if ~exist(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC'), 'file')
            mkdir(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC'))
        end
        save(fullfile(VarCell{CurrentCell}.Experiment.Folder,'HPC',FrameName) , 'VarFrame');
        if F == 1
            Now = clock;
            TimeStamp = [num2str(Now(1)),num2str(Now(2),'%02.0f'),num2str(Now(3),'%02.0f'),'_',num2str(Now(4),'%02.0f'),num2str(Now(5),'%02.0f') ];
            
            VarListFilepath = fullfile(VarCell{CurrentCell}.Experiment.Folder,[TimeStamp, '_VARList.txt']);
            
            fileID = fopen(VarListFilepath, 'a+','n','UTF-8');
            
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
            
            %Write JobArray script file
            JobID = [VarFrame.Experiment.User, '_EN', num2str(VarFrame.Experiment.ExpNum)];
            HPC_ScriptWrite(JobID, TotalNumFrame, fullfile(VarFrame.Analysis.MacPath,[TimeStamp, '_VARList.txt']), VarCell{CurrentCell}.Experiment.Folder)
            
        elseif F == TotalNumFrame
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
            
            fclose(fileID);
            
        else
            fprintf(fileID,'%s\n', fullfile(VarFrame.Analysis.MacPath, 'HPC', FrameName));
        end
        
        
        
    end
elseif ~isempty(strfind(VarRec.Analysis.System, 'Linux'))
    fprintf('linux')
    %Change FilePath
    for C = 1:length(VarCell)
        VarCell{C} = ChangPathForMount(VarCell{C}, LinuxPath);
    end
    %Save VarCell with modified path
    [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
    CellName = ['VarCell_', CellName, '.mat'];
    assignin('base','VarCell',VarCell);
    save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
    
elseif strcmp(VarRec.Analysis.MacServerBase, 'local')
    
    
    %Save VarCell with same path as the one use for ANalysis preparation
    [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
    CellName = ['VarCell_', CellName, '.mat'];
    save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
    
    
else
    %If required modify path names in VarCell to appropriate platform based on
    %ServerBase name
    if strcmpi(System, 'mac') % || strcmpi(System, 'mac_para')
        %For analysis on Mac save VarCell in first experiment folder
        if ismac
            [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
            CellName = ['VarCell_', CellName, '.mat'];
            save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
        else ispc
            %Modify filenames for images to adjust path from PC to Mac
            for C = 1:length(VarCell)
                
                VarCell{C}.Analysis.OutPath = PC2MAC(VarCell{C}.Analysis.OutPath, VarCell{C}.Analysis.MacServerBase);
                
                %Replace path for FilePath variable
                if iscell(VarCell{C}.Analysis.FilePath)
                    [M, N] = size(VarCell{C}.Analysis.FilePath);
                    for m = 1:M
                        for n = 1:N
                            VarCell{C}.Analysis.FilePath{m,n} = PC2MAC(VarCell{C}.Analysis.FilePath{m,n}, VarCell{C}.Analysis.MacServerBase);
                        end
                    end
                else
                    VarCell{C}.Analysis.FilePath = PC2MAC(VarCell{C}.Analysis.FilePath, VarCell{C}.Analysis.MacServerBase);
                end
                
            end
            [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
            CellName = ['VarCell_', CellName, '.mat'];
            save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
        end
        
    elseif strcmpi(System, 'PC') % || strcmpi(System, 'PC_para')
        %For analysis on PC save VarCell in first experiment folder
        if ispc
            [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
            CellName = ['VarCell_', CellName, '.mat'];
            save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
        elseif ismac
            for C = 1:length(VarCell)
                VarCell{C}.Analysis.OutPath = Mac2PC(VarCell{C}.Analysis.OutPath, VarCell{C}.Analysis.PCServerBase);
                
                %Replace path for FilePath variable
                if iscell(VarCell{C}.Analysis.FilePath)
                    [M, N] = size(VarCell{C}.Analysis.FilePath);
                    for m = 1:M
                        for n = 1:N
                            VarCell{C}.Analysis.FilePath{m,n} = Mac2PC(VarCell{C}.Analysis.FilePath{m,n}, VarCell{C}.Analysis.PCServerBase);
                        end
                    end
                else
                    VarCell{C}.Analysis.FilePath = Mac2PC(VarCell{C}.Analysis.FilePath, VarCell{C}.Analysis.PCServerBase);
                end
            end
            [Path, CellName] = fileparts(VarCell{1}.Experiment.Folder);
            if isempty(CellName)
                [Path, CellName] = fileparts(Path);
            end
            assignin('base', 'VarCell', VarCell)
            CellName = ['VarCell_', CellName, '.mat'];
            save(fullfile(VarCell{1}.Experiment.Folder,CellName) , 'VarCell');
        end
    end
end
% catch LastErr
%     VarCell{1}.error = LastErr;
%     assignin('base','VarCell_Error',VarCell);
%     %Save VarRec file
%     Time = clock;
%     TimeStr = [num2str(Time(1)),num2str(Time(2),'%02d'),num2str(Time(3),'%02d'),'_',num2str(Time(4),'%02d'),num2str(Time(5),'%02d')];
%     save(fullfile(BaseDir,['ERROR_',TimeStr,'.mat']), 'VarCell')
% end





%% Transform MAC to PC files %%%%%%%%%%%%%%%%%%%

function PCOut = Mac2PC(Path, PCbase)

%Replace FilePath and OutPath with PC server path
%Find PC Folder name (Should match one of the folder name in the
%original path
if ~isempty(Path)
    if iscell(Path)
        Sep = strfind(PCbase, '\');
        NewRoot = PCbase(1:Sep(end));
        if Sep(end) == length(PCbase)
            CommonFolder = PCbase(Sep(end-1)+1:end-1);
        else
            CommonFolder = PCbase(Sep(end)+1:end);
        end
        NumZ = length(Path);
        for Z = 1:NumZ
            
            %Go up the path to find the same folder in the Mac path
            [OldRoot, UpFolder] = fileparts(Path{Z});
            %Loop through all folders until match is found between PC and mac
            %folder name
            while ~strcmp(UpFolder, CommonFolder)
                PrevRoot = OldRoot;
                [OldRoot, UpFolder] = fileparts(PrevRoot);
                if isempty(UpFolder)
                    error(['Error generating MAC path: No Match found between ',Path{Z} 'and  ', PCbase])
                end
            end
            %Replace path for OutPath variable
            PCOut{Z} = [NewRoot,Path{Z}(length(PrevRoot)+2:end)];
            PCOut{Z} = regexprep(PCOut{Z}, '/', '\');
        end
        
    else
        
        Sep = strfind(PCbase, '\');
        NewRoot = PCbase(1:Sep(end));
        if Sep(end) == length(PCbase)
            CommonFolder = PCbase(Sep(end-1)+1:end-1);
        else
            CommonFolder = PCbase(Sep(end)+1:end);
        end
        %Go up the path to find the same folder in the Mac path
        [OldRoot, UpFolder] = fileparts(Path);
        %Loop through all folders until match is found between PC and mac
        %folder name
        while ~strcmp(UpFolder, CommonFolder)
            PrevRoot = OldRoot;
            [OldRoot, UpFolder] = fileparts(PrevRoot);
            if isempty(UpFolder)
                error(['Error generating MAC path: No Match found between ',Path 'and  ', PCbase])
            end
        end
        %Replace path for OutPath variable
        PCOut = [NewRoot,Path(length(PrevRoot)+2:end)];
        PCOut = regexprep(PCOut, '/', '\');
    end
else
    PCOut = [];
    
end



%% Transform PC to MAC files %%%%%%%%%%%%%%%%%%%



function MacOut = PC2MAC(Path, Macbase)

%Replace FilePath and OutPath with PC server path
%Find PC Folder name (Should match one of the folder name in the
%original path
if ~isempty(Path)
    
    Sep = strfind(Macbase, '/');
    NewRoot = Macbase(1:Sep(end));
    
    if Sep(end) == length(Macbase)
        CommonFolder = Macbase(Sep(end-1)+1:end-1);
    else
        CommonFolder = Macbase(Sep(end)+1:end);
    end
    %Go up the path to find the same folder in the Mac path
    [OldRoot, UpFolder] = fileparts(Path);
    %Loop through all folders until match is found between PC and mac
    %folder name
    while ~strcmp(UpFolder, CommonFolder)
        PrevRoot = OldRoot;
        [OldRoot, UpFolder] = fileparts(PrevRoot);
        if isempty(UpFolder)
            error(['Error generating PC path: No Match found between ',Path 'and  ', Macbase])
        end
    end
    %Replace path for OutPath variable
    MacOut = [NewRoot,Path(length(PrevRoot)+2:end)];
    MacOut = regexprep(MacOut, '\', '/');
else
    MacOut = [];
    
end


function VarOUT = ChangPathForWally(VarIN, NewPath)

%OldPath = fileparts(fileparts(VarIN.Analysis.MacPath));
StrInd = strfind(VarIN.Analysis.MacPath, 'spelet1');
OldPath = VarIN.Analysis.MacPath(1:StrInd+length('spelet1')-1);
%Duplicate VarIN
VarOUT = VarIN;

%Replace Back slashes:
NewPath = regexptranslate('escape', NewPath);
OldPath = regexptranslate('escape', OldPath);

%loop through all file names to modify them
for i = 1:size(VarIN.Analysis.FilePath,1)
    for t = 1:size(VarIN.Analysis.FilePath,2)
        %         NumZ =  VarIN.Analysis.Zstack
        if isfield(VarIN.Analysis, 'Zstack') && VarIN.Analysis.Zstack(i) >1 && ~isempty(VarIN.Analysis.FilePath{i,t})
            for z = 1:VarIN.Analysis.Zstack(i)
                %Get old file
                OldFile = VarIN.Analysis.FilePath{i,t}{z};
                %Change directory
                NewFile = regexprep(OldFile, OldPath, NewPath);
                %Add file to structure
                VarOUT.Analysis.FilePath{i,t}{z} = NewFile;
            end
        elseif ~isempty(VarIN.Analysis.FilePath{i,t})
            
            %Get old file
            OldFile = VarIN.Analysis.FilePath{i,t};
            
            %Change directory
            NewFile = regexprep(OldFile, OldPath, NewPath);
            
            %Change file separator
            %   if strcmpi(NewSystem, 'Mac')
            %             FS = regexptranslate('escape', '\');
            %        NewFile = regexprep(NewFile, FS, '/');
            %Add file to structure
            VarOUT.Analysis.FilePath{i,t} = NewFile;
        else
            VarOUT.Analysis.FilePath{i,t} = [];
        end
        
        if i == 1&& t == 1
            VarOUT.Analysis.MacPath = regexprep(VarIN.Analysis.MacPath, OldPath, NewPath);
            %VarOUT.Analysis.MacPath = regexprep( VarOUT.Analysis.MacPath, FS, '/');
            
            VarOUT.Analysis.OutPath = regexprep(VarIN.Analysis.OutPath, OldPath, NewPath);
            if isfield(VarOUT.Analysis, 'FlatFile')
                VarOUT.Analysis.FlatFile = regexprep(VarIN.Analysis.FlatFile, OldPath, NewPath);
            end
            
        end
    end
end



function VarOUT = ChangPathForCluster(VarIN, NewPath)

%OldPath = fileparts(fileparts(VarIN.Analysis.MacPath));
StrInd = strfind(VarIN.Analysis.MacPath, 'gr_Pelet');
OldPath = VarIN.Analysis.MacPath(1:StrInd+length('gr_Pelet')-1);
%Duplicate VarIN
VarOUT = VarIN;

%Replace Back slashes:
NewPath = regexptranslate('escape', NewPath);
OldPath = regexptranslate('escape', OldPath);

%loop through all file names to modify them
for i = 1:size(VarIN.Analysis.FilePath,1)
    for t = 1:size(VarIN.Analysis.FilePath,2)
        %         NumZ =  VarIN.Analysis.Zstack
        if isfield(VarIN.Analysis, 'Zstack') && VarIN.Analysis.Zstack(i) >1 && ~isempty(VarIN.Analysis.FilePath{i,t})
            for z = 1:VarIN.Analysis.Zstack(i)
                %Get old file
                OldFile = VarIN.Analysis.FilePath{i,t}{z};
                %Change directory
                NewFile = regexprep(OldFile, OldPath, NewPath);
                %Add file to structure
                VarOUT.Analysis.FilePath{i,t}{z} = NewFile;
            end
        elseif ~isempty(VarIN.Analysis.FilePath{i,t})
            
            %Get old file
            OldFile = VarIN.Analysis.FilePath{i,t};
            
            %Change directory
            NewFile = regexprep(OldFile, OldPath, NewPath);
            
            %Change file separator
            %   if strcmpi(NewSystem, 'Mac')
            %             FS = regexptranslate('escape', '\');
            %        NewFile = regexprep(NewFile, FS, '/');
            %Add file to structure
            VarOUT.Analysis.FilePath{i,t} = NewFile;
        else
            VarOUT.Analysis.FilePath{i,t} = [];
        end
        
        if i == 1&& t == 1
            VarOUT.Analysis.MacPath = regexprep(VarIN.Analysis.MacPath, OldPath, NewPath);
            %VarOUT.Analysis.MacPath = regexprep( VarOUT.Analysis.MacPath, FS, '/');
            
            VarOUT.Analysis.OutPath = regexprep(VarIN.Analysis.OutPath, OldPath, NewPath);
            if isfield(VarOUT.Analysis, 'FlatFile')
                VarOUT.Analysis.FlatFile = regexprep(VarIN.Analysis.FlatFile, OldPath, NewPath);
            end
            
        end
    end
end


function VarOUT = ChangPathForMount(VarIN, NewRoot)

OldRoot = VarIN.Analysis.MacServerBase;
%Duplicate VarIN
VarOUT = VarIN;

%Replace Back slashes:
NewRoot = regexptranslate('escape', NewRoot);
OldRoot = regexptranslate('escape', OldRoot);

%loop through all file names to modify them
for i = 1:size(VarIN.Analysis.FilePath,1)
    for t = 1:size(VarIN.Analysis.FilePath,2)
        %         NumZ =  VarIN.Analysis.Zstack
        if isfield(VarIN.Analysis, 'Zstack') && VarIN.Analysis.Zstack(i) >1 && ~isempty(VarIN.Analysis.FilePath{i,t})
            for z = 1:VarIN.Analysis.Zstack(i)
                %Get old file
                OldFile = VarIN.Analysis.FilePath{i,t}{z};
                %Change directory
                NewFile = regexprep(OldFile, OldRoot, NewRoot);
                %Add file to structure
                VarOUT.Analysis.FilePath{i,t}{z} = NewFile;
            end
        elseif ~isempty(VarIN.Analysis.FilePath{i,t})
            
            %Get old file
            OldFile = VarIN.Analysis.FilePath{i,t};
            
            %Change directory
            NewFile = regexprep(OldFile, OldRoot, NewRoot);
            
            %Change file separator
            %   if strcmpi(NewSystem, 'Mac')
            %             FS = regexptranslate('escape', '\');
            %        NewFile = regexprep(NewFile, FS, '/');
            %Add file to structure
            VarOUT.Analysis.FilePath{i,t} = NewFile;
        else
            VarOUT.Analysis.FilePath{i,t} = [];
        end
        
        if i == 1&& t == 1
            VarOUT.Analysis.MacPath = regexprep(VarIN.Analysis.MacPath, OldRoot, NewRoot);
            %VarOUT.Analysis.MacPath = regexprep( VarOUT.Analysis.MacPath, FS, '/');
            
            VarOUT.Analysis.OutPath = regexprep(VarIN.Analysis.OutPath, OldRoot, NewRoot);
            if isfield(VarOUT.Analysis, 'FlatFile')
                VarOUT.Analysis.FlatFile = regexprep(VarIN.Analysis.FlatFile, OldRoot, NewRoot);
            end
            
        end
    end
end

