function ParallelFrameAnalysis(VarCell, OverWrite)
%VarCell can either be the cell array to analyze or a string representing
%the path where the data is loaded

%FOR CellX analysis: Use Path to VarCell. Otherwise CellXObject is not
%recognized when loading the VarCell in Matlab without the proper paths

%Overwrite = -1: Prepare only the VarFrame data
%OverWrite = 0: Perform segmentation only on missing Frame.mat files 
%OverWrite = 1: Delete Data.mat file and reprocess the export file
%OverWrite = 2: Delete Frame.mat file and re-do the segmentation

BaseDir = cd;
addpath(fullfile(BaseDir, 'CommonProcess'))
addpath(genpath('CellX/core'));

close all

if nargin == 1
    OverWrite = 0;
end

%Check if VarCell is a cell
if ~iscell(VarCell)
    VarCellPath = VarCell;
    load(VarCellPath)
end

%Calculate the total number of frames to be analyzed
NumCells = length(VarCell);
TotalNumFrame = 0;
for P = 1:NumCells
    if isfield(VarCell{P}.Experiment, 'Skip')
        VarCell{P}.Analysis.NumSkip = max(VarCell{P}.Analysis.SkipFrame(1,:));
        TotalNumFrame = TotalNumFrame + ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
        VarCell{P}.Analysis.NumFrame2Analyze = ceil(VarCell{P}.Analysis.NumFrame/VarCell{P}.Analysis.NumSkip);
    else
        TotalNumFrame = TotalNumFrame + VarCell{P}.Analysis.NumFrame;
        VarCell{P}.Analysis.NumFrame2Analyze = VarCell{P}.Analysis.NumFrame;
        VarCell{P}.Analysis.NumSkip = 1;
    end
    
%      OutPath = strrep(VarCell{P}.Analysis.OutPath, '/mnt/vwosika/', '/Volumes/DMF/GROUPS/gr_Pelet/')
%     VarCell{P}.Analysis.OutPath = OutPath;
    
    ImgFolder = fileparts(VarCell{P}.Analysis.OutPath);
    MATOutFolder = fullfile(ImgFolder,['MATOut_', num2str(VarCell{P}.Analysis.CurrentPos, '%03d')]);
    if exist(MATOutFolder, 'file') == 0
        mkdir(MATOutFolder)
    elseif OverWrite == 2
        %fprintf('delete MATOut')
        rmdir(MATOutFolder,'s')
        mkdir(MATOutFolder)
    end
    DataOutFolder = fullfile(ImgFolder,'DATAOut');
    if exist(DataOutFolder, 'file') == 0
        mkdir(DataOutFolder)
    elseif OverWrite >= 1
        %fprintf('delete DataOut')
        rmdir(DataOutFolder,'s')
        mkdir(DataOutFolder)
    end
    
    %Generate Export filename
    [SaveDir] = fileparts(ImgFolder);
    [~, ExperimentFolderName] = fileparts(SaveDir);
    
    ExportFileName = fullfile(SaveDir,['Export_', ExperimentFolderName, '.mat']);
    %Delete export file
    if exist(ExportFileName, 'file')
        delete(ExportFileName)
    end
end




 % Calculate Position and Frame number to generate VarFrame Cell array
for F = 1:TotalNumFrame
    CurrentFrame = F;
    CurrentCell = 1;
    while CurrentCell<NumCells && CurrentFrame > VarCell{CurrentCell}.Analysis.NumFrame2Analyze
        CurrentFrame = CurrentFrame - VarCell{CurrentCell}.Analysis.NumFrame2Analyze;
        CurrentCell = CurrentCell+1;
    end
    
     %If Skipped frames:
    if VarCell{CurrentCell}.Analysis.NumSkip >1
        CurrentFrame = VarCell{CurrentCell}.Analysis.NumSkip*(CurrentFrame-1)+1;
    end
    %Transfer Position and Frame numbers to Var
    VarFrame{F} = VarCell{CurrentCell};
    VarFrame{F}.Analysis.CurrentCell = CurrentCell;
    VarFrame{F}.Analysis.CurrentFrame = CurrentFrame;

    %Check if Skipped frame for analysis
    if CurrentFrame > 1 && isfield(VarFrame{F}.Experiment, 'Skip') && min(VarFrame{F}.Analysis.SkipFrame(1,:))< 0
        
        SkipAnalysis = find(VarFrame{F}.Analysis.SkipFrame(1,:)<0);
        for SA = 1:length(SkipAnalysis)
            VarFrame{F}.Analysis.SkipFrame(2,SkipAnalysis(SA)) =  VarFrame{F}.Analysis.SkipFrame(2,SkipAnalysis(SA)) + 1;
            if VarFrame{F}.Analysis.SkipFrame(2,SkipAnalysis(SA)) > abs(VarFrame{F}.Analysis.SkipFrame(1,SkipAnalysis(SA)))
                VarFrame{F}.Analysis.SkipFrame(2,SkipAnalysis(SA)) = 1;
            end
             VarCell{CurrentCell}.Analysis.SkipFrame(2,SkipAnalysis(SA)) = VarFrame{F}.Analysis.SkipFrame(2,SkipAnalysis(SA));
        end
    end
    
end
assignin('base', 'VarFrame',VarFrame)
if OverWrite >= 0
    %% Start the workers if running in parallel mode
    NumWorkers = 8;
    
    %Set number of workers to be activated
    if TotalNumFrame>= NumWorkers
        parpool(NumWorkers)
    else
        parpool(TotalNumFrame)
    end
    
    
    
    %% Perform Frame analysis
    parfor F = 1:TotalNumFrame
        %for F = 1:TotalNumFrame
        
        FrameAnalysis(VarFrame{F})
    end
    
end