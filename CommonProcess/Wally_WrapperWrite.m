function Wally_WrapperWrite(JobID, TimeStamp, VarCell)

% Writte Wrapper script for submission and Data retrieval

%Create Wrapper File for submission
NASBasefolderPath = VarCell{1}.Experiment.Folder;
UserName = VarCell{1}.Analysis.RunningUser;

WrapperPath = fullfile(NASBasefolderPath, ['Wrapper_',JobID, '.sh']);
fileID= fopen(WrapperPath, 'w+');


for C = 1: numel(VarCell)
     SplitPath = split(VarCell{C}.Experiment.Folder, filesep);
     FACPath = strjoin(SplitPath(4:end), filesep);
    
   
     
     if C == 1
          %Create directories
        MkDir_Cmd = ['mkdir -p /scratch/wally/', FACPath];
        fprintf(fileID,'%s\n',MkDir_Cmd);
          
           %MoveData
         
            SyncDat_Cmd = ['syncdat -r /scratch/wally/', FACPath, ' ',UserName, '@dexwin.unil.ch:' FACPath];
            fprintf(fileID,'%s\n',SyncDat_Cmd);
            PrevFacPath = FACPath;
            Iter = 1;
            AllFACPath{Iter} = FACPath;
            
            %Crate output directory
          MkDir_Cmd = ['mkdir -p /scratch/wally/', FACPath '/', TimeStamp, '_OUT'];
          fprintf(fileID,'%s\n',MkDir_Cmd);
     else
         if ~strcmp(PrevFacPath,FACPath)
              %Create directories
        MkDir_Cmd = ['mkdir -p /scratch/wally/', FACPath];
        fprintf(fileID,'%s\n',MkDir_Cmd);
        
        %MoveData
            SyncDat_Cmd = ['syncdat -r /scratch/wally/', FACPath, ' ',UserName, '@dexwin.unil.ch:' FACPath];
            fprintf(fileID,'%s\n',SyncDat_Cmd);
            PrevFacPath = FACPath;
             Iter = Iter+1;
            AllFACPath{Iter} = FACPath;
         end
     end
     
     

    
end

%sbatch command

sbatch_Cmd = ['sbatch /scratch/wally/', FACPath, '/batchArray_', JobID,'.sh'];
fprintf(fileID,'%s\n',sbatch_Cmd);

fclose(fileID);


%% Wrapper for Export and DataOUT retrieval

%Create Wrapper File 


WrapperPath = fullfile(NASBasefolderPath, ['SendBack_',JobID, '.sh']);
fileID= fopen(WrapperPath, 'w+');


for P = 1: numel(AllFACPath)
    Move_Cmd = ['movedat -i -B /scratch/wally/', AllFACPath{P}, '/DATAOut ',UserName, '@dexwin.unil.ch:' , AllFACPath{P}, '/DATAOut'];
     fprintf(fileID,'%s\n',Move_Cmd);
    Move_Cmd = ['movedat /scratch/wally/', AllFACPath{P}, '/Export_*.mat ',UserName, '@dexwin.unil.ch:' , AllFACPath{P}, '/'];
    fprintf(fileID,'%s\n',Move_Cmd);
end
fclose(fileID);

