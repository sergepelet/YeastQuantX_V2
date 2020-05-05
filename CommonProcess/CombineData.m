function Var = CombineData(Var)

tic
%Get Number of position
NumPos = Var.Analysis.NumPos;

%Get Saving directory
[SaveDir] = fileparts(Var.Analysis.OutPath);
[~, FolderName] = fileparts(SaveDir);
% SaveDir = strrep(SaveDir, '/mnt/spelet1/', '/Volumes/DMF/GROUPS/gr_Pelet/')
%SaveDir = '/Users/serge/Projects/ImageAnalysis/MeasAnalysis/2016/160428_DA';
ExpFileName = ['Export_', FolderName, '.mat'];
% Test if export file exists

fprintf('\nCreate Export File\n')
Export = [];
Var.Analysis.ExportFileCreated = fullfile(SaveDir, ExpFileName);


NotExported = [];
for P = 1:NumPos
    ImgFolder = fileparts(Var.Analysis.OutPath);
   %test if Data.mat file can be loaded
    try
        DataFileName = fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(P, '%03d'),'_Data.mat']);
        load(DataFileName);
    catch LoadError
        fprintf(['Pos_',num2str(P, '%03d'),'_Data.mat not loaded in export\n'])
        NotExported = [NotExported, P];
        continue
    end
    if ~isfield(Var, 'error')
        %Get name of condition
        CondName = Var.Analysis.CurrentCondition;  %['Well_',num2str(Var.Analysis.CurrentWell)];%
        %Corect for non authorized character in string name
        CondName = regexprep(CondName, ' ', '_');
        CondName= regexprep(CondName, '\.', 'p');
        CondName= regexprep(CondName, '-', '_');
        CondName = regexprep(CondName, '\s', '_');
        
        if ~isvarname(CondName)
            CondName = ['c_',CondName];
        end
        
        %Generate time axis
        TimeAxis =Var.Analysis.TimeAxis./1000;
        %Make sure that time increases linearly
        %for Experiments that take place over midnight Add 24 hours shift to
        %time stamp
        for TA = 2: length(TimeAxis)
            if TimeAxis(TA) < TimeAxis(TA-1)
                TimeAxis(TA) = TimeAxis(TA)+24*60*60;
            end
        end
        %Get TimeStamp
        TimeStamp =Var.Analysis.TimeStamp';
        %Get names of all objects
        ObjFields = fieldnames(Var.Measurements);
        %Get names of all illuminations
        IllumNames = Var.Analysis.MeasureIntensity;
        
        
        %Check if only track cells are exported
        if isfield(Var.Analysis, 'ExportType') && (strcmp(Var.Analysis.ExportType, 'Tracked'))
            %Get indices of tracked cells
            ExportedCell = find(min(Var.Measurements.(ObjFields{1}).(IllumNames{1}).CheckLabel,[],2)>0);
        else
            %Get indices of all cells
            ExportedCell = 1:Var.Analysis.NumberCell;
        end
        
        %Loop through all Objects
        for Obj = 1:length(ObjFields)
            %Get Measurement list
            if isfield(Var.Analysis, 'ExportAllMeas') && strcmp(Var.Analysis.ExportAllMeas, 'Yes')
                MeasList = fieldnames(Var.Measurements.(ObjFields{Obj}).(IllumNames{1}));
                MeasList{length(MeasList)+1} = 'Position';
            else
                MeasList = Var.Analysis.ExportMeasure;
                MeasList{length(MeasList)+1} = 'CheckLabel';
                MeasList{length(MeasList)+1} = 'Position';
                MeasList{length(MeasList)+1} = 'CenterX';
                MeasList{length(MeasList)+1} = 'CenterY';
                MeasList{length(MeasList)+1} = 'Eccentricity';
                MeasList{length(MeasList)+1} = 'MajorAxisLength';
                MeasList{length(MeasList)+1} = 'MinorAxisLength';
                MeasList{length(MeasList)+1} = 'ConnectedHiPix';
                MeasList{length(MeasList)+1} = 'ConnectedHiPixList';
                if isfield(Var.Measurements.(ObjFields{Obj}).(IllumNames{1}), 'ParentLabel')
                    MeasList{length(MeasList)+1} = 'ParentLabel';
                end
            end
            
            %loop through Illum
            for I = 1:length(IllumNames)
                %Check if Export structure already exists for this object
                if isfield(Export, CondName) && isfield(Export.(CondName),IllumNames{I}) ...
                        && isfield(Export.(CondName).(IllumNames{I}),ObjFields{Obj})...
                        && isfield(Export.(CondName).(IllumNames{I}).(ObjFields{Obj}),(MeasList{1}))
                    %loop through all measurements
                    for M = 1:length(MeasList)
                        %Check if Meas is position
                        if strcmp(MeasList{M}, 'Position')
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                [Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells, ...
                                P.*ones(Var.Analysis.NumFrame, length(ExportedCell))];
                        elseif strcmp(MeasList{M}, 'TimeStamp')
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                [Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells, ...
                                TimeStamp*ones(1, length(ExportedCell))];
                        else
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                [Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells, ...
                                Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasList{M})(:,ExportedCell)];
                        end
                    end
                else  %New export structure
                    %loop through all measurements
                    for M = 1:length(MeasList)
                        %Check if Meas is position
                        if strcmp(MeasList{M}, 'Position')
                            if Obj == 1 && I == 1
                                %Transfer Experiment description to export stucture
                                Export.(CondName).Experiment = Var.Experiment;
                                %Transfer time axis to Export
                                Export.(CondName).Time = TimeAxis;
                            end
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                P.*ones(Var.Analysis.NumFrame, length(ExportedCell));
                        elseif strcmp(MeasList{M}, 'TimeStamp')
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                TimeStamp*ones(1, length(ExportedCell));
                            
                        else
                            Export.(CondName).(IllumNames{I}).(ObjFields{Obj}).(MeasList{M}).Cells = ...
                                Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasList{M})(:,ExportedCell);
                        end
                        
                        
                    end
                end
            end
        end
        
    end
end
if ~isempty(NotExported)
    Export.MissingPositions.PositionNumber = NotExported;
    Export.MissingPositions.LastLoadError = LoadError;
end
assignin('base','Export',Export);
%Save export variable
save(fullfile(SaveDir, ExpFileName), 'Export', '-v7.3')
%Save Timing Info
Var.Analysis.Timing.(mfilename) = toc;