function FrameAnalysis(Var)

close all

%Check if VAR Frame_XXX.mat exists
for S = 1:Var.Analysis.NumSkip
    ImgFolder = fileparts(Var.Analysis.OutPath);
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame+S-1,'%05d'),'.mat']);
    NoMATList(S) = exist(FileName, 'file');
end
NoMAT = min(NoMATList);

%If file exists Verify that it does not contain an error field
if NoMAT~= 0
    for S = 1:Var.Analysis.NumSkip
        ImgFolder = fileparts(Var.Analysis.OutPath);
        FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame+S-1,'%05d'),'.mat']);
        SavedVar = load(FileName);
        if isfield(SavedVar.Var, 'error')
            delete(FileName)
            NoMAT = 0;
            clear SavedVar
        end
    end
end


    

if NoMAT == 0
 %   assignin('base','Var_start',Var);
    %% Individual Frame loading and segmentation
  try
        %Perfrom segmentation on first frame
        for S = 1:Var.Analysis.NumSkip
      %     S = S
      %    SF = Var.Analysis.SkipFrame
            if S >1
             %   assignin('base',['Var_start', num2str(S)],VarNext);
                %Remove Image measurements from Structure
                Var = RemoveMeasurements(VarNext);
            end
            %% Load Images
            
            for L = 1:length(Var.Analysis.LoadIllum)
                Var = LoadImage(Var, L);
                drawnow
                
                %Var = CorrectFlatness(Var, L);
            end
           
            if S == 1
                Var.Analysis.SegmentedFrame = Var.Analysis.CurrentFrame;
                %% Segmentation routine
                %assignin('base','Var_BeforeSeg',Var);
                if strcmp(Var.Experiment.Segmentation, 'Nucl_BF')
                    %Segment Nucleus in images
                    %Var = SegmentYeastFluo(Var, 1);
                     Var = SegmentYeastFluo(Var, 1);
                   % assignin('base',num2str(Var.Analysis.CurrentFrame,'Var_SegFluo_%03d'),Var);
                    %Get cell area around segmented Nuclei with Brightfield image
                    assignin('base','Var_BeforeBF',Var);
                    Var = SegmentNuclBF(Var, 1);
                elseif strcmp(Var.Experiment.Segmentation, 'CellXNucl')
                    %Segment Nucleus in images
                     Var = SegmentYeastFluo(Var, 1);
                     %Get cell area around segmented Nuclei with Brightfield image
                  assignin('base','Var_BeforeBF',Var);
                    Var = SegmentCellX(Var, 1);
                   % error('bla')
                elseif strcmp(Var.Experiment.Segmentation, 'FluoMedia')  %% TO BE TESTED
                    %Segment cell based on fluorescence from media in microfluidic chips
                    Var = SegmentYeastFluoMedia(Var, 1);
                    
                elseif strcmp(Var.Experiment.Segmentation, 'FluoCell') %% TO BE TESTED
                    %Segment Nucleus or cell based on Fluorescence intensity
                    Var = SegmentYeastFluo(Var, 1);
                elseif strcmp(Var.Experiment.Segmentation, 'BrightField') %% TO BE TESTED
                    
                    %Segment DIC images
                    Var = SegmentYeastPhase(Var, 1);
                    %Find individual cells using Hough transform
                    Var = SplitGroupHough(Var, 1);
                elseif strcmp(Var.Experiment.Segmentation, 'mRNA')
                    %Segment Nucleus in images
                    Var = SegmentYeastFluo(Var, 1);
                    
                    %Get cell area around segmented Nuclei with Brightfield image
                    Var = SegmentNuclBF(Var, 1);
                    assignin('base','Var_NBF',Var);
                    
                    %Segment Dots in cells
                    %Var = SegmentDots(Var,1);
                end
              %assignin('base','Var_Seg',Var);
                % error('bla')
              
            else
%                 assignin('base','Var_BOT',Var);
%                 Var = ObjectTransfer(Var);
                
                if strcmp(Var.Experiment.Segmentation, 'mRNA')
                    %Segment Dots in cells
                    Var = SegmentDots(Var,1);
                    
                end
        
            end
            
            %% Secondary Object definition
            if isfield(Var.Analysis, 'SecPrimObj')
                for L = 1:length(Var.Analysis.SecPrimObj)
                    Var = SecondaryObject(Var, L);
                end
            end
            
          % assignin('base','Var_2',Var);
            if isfield(Var.Analysis, 'ExpObjSmall')
                for L = 1:length(Var.Analysis.ExpObjSmall)
                    Var = ExpandObject(Var, L);
                end
            end
        %    assignin('base','Var_SecObj',Var);
            %% Object measurements
            
            %Subtract background from Intensity Image
            if isfield(Var.Analysis, 'SubtractIn')
                for L = 1:length(Var.Analysis.SubtractIn)
                    Var = SubtractBackground(Var, L);
                end
            end
      assignin('base','Var_BM',Var);
            %Measure Objects
            if isfield(Var.Analysis, 'MeasureIntensity') && isfield(Var.Analysis, 'MeasureObject')
                for I = 1:length(Var.Analysis.MeasureIntensity)
                    for O = 1:length(Var.Analysis.MeasureObject)
                        L = I+O/100;
                        Var = MeasureObject(Var, L);
                    end
                end
            end
        %    assignin('base','Var_Meas',Var);
            %Calculate new object from measurements
            if isfield(Var.Analysis, 'CalcObjOUT')
                for L = 1:length(Var.Analysis.CalcObjOUT)
                    Var = CalculatedObject(Var, L);
                end
            end
            
            % Make movie
            if isfield(Var.Analysis, 'MovieImg')
                for L = 1:length(Var.Analysis.MovieImg)
                    Var = Img2MovieRGB(Var, L);
                end
            end
            
            if S < Var.Analysis.NumSkip
                VarNext = Var;
            end
            
            %% Save VAR to frame.mat file
            %Remove Img field from Data
            Var = rmfield(Var, 'Img');
            %Save measurements to disk.
            ImgFolder = fileparts(Var.Analysis.OutPath);
            FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame,'%05d'),'.mat']);
            save(FileName, 'Var', '-v7.3');
                
            %Generate Verif TXT file
            TxtFrame = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame,'%05d'),'.txt']);
            
            fileID = fopen(TxtFrame,'w');
            fprintf(fileID, ['Frame analysis', num2str(Var.Analysis.CurrentFrame,'%05d DONE')]);
            fclose(fileID);
            
            if S < Var.Analysis.NumSkip
                %Increment Current Frame by one
                VarNext.Analysis.CurrentFrame = VarNext.Analysis.CurrentFrame+1;
                %% Current Frame is larger than NumFrame then break the loop
                if VarNext.Analysis.CurrentFrame > VarNext.Analysis.NumFrame
                    break
                end
            end
            
            
            
        end
    catch LastErr
        
        
        if isfield(Var.Analysis, 'CurrentFileName')
            fprintf(['\nERROR IN ', Var.Analysis.CurrentFileName, ' Pos: ', num2str(Var.Analysis.CurrentPos),...
                ' Frame: ', num2str(Var.Analysis.CurrentFrame),'\n'])
        else
            fprintf('\nERROR in ImgLoad: Check the image filepath\n')
        end
        
        
        ErrorString = [];
        for E = 1:length(LastErr.stack)
            ErrorString = [ErrorString, '\nIn function ', LastErr.stack(E,1).name, ' at line ', num2str(LastErr.stack(E,1).line)];
        end
        fprintf(ErrorString)
        fprintf('\n')
        LastErrMessage = LastErr.message;
        fprintf(LastErrMessage)
        fprintf('\n')
        fprintf('\n')
        Var.error = LastErr;
        
        assignin('base','Var_Error',Var);
        %Save Frame.mat with error field
        ImgFolder = fileparts(Var.Analysis.OutPath);
        FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame,'%05d'),'.mat']);
        save(FileName, 'Var');
        
         ErrorText = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Error_F', num2str(Var.Analysis.CurrentFrame,'%05d'),'.txt']);
        fileID = fopen(ErrorText,'w');
        fprintf(fileID, ErrorString);
        fprintf(fileID,'\n');
        fprintf(fileID,LastErrMessage);
        fclose(fileID);
        
        %Generate Verif TXT file
            TxtFrame = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.CurrentFrame,'%05d'),'.txt']);
            
            fileID = fopen(TxtFrame,'w');
            fprintf(fileID, ['Frame analysis', num2str(Var.Analysis.CurrentFrame,'%05d Finished with error')]);
            fclose(fileID);
        
        %Save Data Mat as only the current Frame.mat Var with error field
        DataFileName = fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(Var.Analysis.CurrentPos, '%03d'),'_Data.mat']);
        
        save(DataFileName, 'Var', '-v7.3');
        
         %Generate Verif TXT file for Data.mat
            TxtData =  fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(Var.Analysis.CurrentPos, '%03d'),'_Data.txt']);
            
            fileID = fopen(TxtData,'w');
            fprintf(fileID, ['Data.mat for Pos ', num2str(Var.Analysis.CurrentPos,'%05d with error')]);
            fclose(fileID);
        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Process frame.mat files              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get list of MAT Frame files
MATfolder = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')]);
TxtFrameList = dir(fullfile(MATfolder, 'Frame_*.txt'));
NumFrameFiles = length(TxtFrameList);

%Create DataFile name
ImgFolder = fileparts(Var.Analysis.OutPath);
DataFileName = fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(Var.Analysis.CurrentPos, '%03d'),'_Data.mat']);


% Check if Number of MAT file is equal to number of frame
%If it is the case all frames have been analyzed => Can further process data
%Also Check if Data.mat file alread exist
try
    if Var.Analysis.NumFrame == length(TxtFrameList) && exist(DataFileName, 'file') == 0
        
        %Load Frame.mat data
        Var = LoadFrameMeas(Var);
        assignin('base','Var_LoadFr',Var);
        %Track Main object
        Var = TrackObjectsSingle(Var, 1);
        
        assignin('base','VarTrack',Var);
        
        %Combine Frame measurements in single matrix
        Var = CombineFrame(Var);
        
        assignin('base','VarComb',Var);
        % Save Var to Data.mat file
        save(DataFileName, 'Var');
        
    
        %Generate Verif TXT file for Data.mat
            TxtData =  fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(Var.Analysis.CurrentPos, '%03d'),'_Data.txt']);
            
            fileID = fopen(TxtData,'w');
            fprintf(fileID, ['Data.mat for Pos ', num2str(Var.Analysis.CurrentPos,'%05d DONE')]);
            fclose(fileID);
        
        %% Combine Data.mat to create Export
        
        %Check Number of Data. mat files created
        TxtDataList = dir(fullfile(ImgFolder,'DATAOut', 'Pos_*.txt'));
        NumDataFiles = length(TxtDataList);
        %If same number of Data.mat file as position: create Export
        if Var.Analysis.NumPosAnalyzed == NumDataFiles
            fprintf('CombineData\n')
            %Generate Export file
            Var = CombineData(Var);
        end
    else
        %Check Number of Data.mat files
        TxtDataList = dir(fullfile(ImgFolder,'DATAOut', 'Pos_*.txt'));
        NumDataFiles = length(TxtDataList);
        %Generate Export filename
        [SaveDir] = fileparts(ImgFolder);
        [~, FolderName] = fileparts(SaveDir);

        ExportFileName = fullfile(SaveDir,['Export_', FolderName, '.mat']);
        
        if Var.Analysis.NumPosAnalyzed == NumDataFiles && exist(ExportFileName, 'file') == 0
            fprintf('CombineData\n')
            %Generate Export file
            Var = CombineData(Var);
        end
    end
catch LastErr
    
     if isfield(Var.Analysis, 'CurrentFileName')
        fprintf(['\nERROR IN ', Var.Analysis.CurrentFileName, ' combination\n'])
     else
        fprintf(['\nERROR IN file combination\n'])
     end
    
    
    
    ErrorString = [];
    for E = 1:length(LastErr.stack)
        ErrorString = [ErrorString, '\nIn function ', LastErr.stack(E,1).name, ' at line ', num2str(LastErr.stack(E,1).line)];
    end
    fprintf(ErrorString)
    fprintf('\n')
    LastErrMessage = LastErr.message;
    fprintf(LastErrMessage)
    fprintf('\n')
    fprintf('\n')
    Var.error = LastErr;
    
    assignin('base','Var_Error',Var);
    %Save Frame.mat with error field
    ImgFolder = fileparts(Var.Analysis.OutPath);
    %Save again Data Mat as only the current Frame.mat Var with error field
    DataFileName = fullfile(ImgFolder,'DATAOut', ['Pos_',num2str(Var.Analysis.CurrentPos, '%03d'),'_Data.mat']);
    
    save(DataFileName, 'Var');
    
        ErrorText = fullfile(ImgFolder,'DATAOut',['Error_P', num2str(Var.Analysis.CurrentPos,'%03d'),'.txt']);
        fileID = fopen(ErrorText,'w');
        fprintf(fileID, ErrorString);
        fprintf(fileID,'\n')
        fprintf(fileID,LastErrMessage)
        fclose(fileID);
end
