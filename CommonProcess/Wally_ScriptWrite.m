function Wally_ScriptWrite(JobID, NumFrame, VarListPath, NASfolderPath)
%Get VARLIST file info


%% Open template for batch array script

BaseDir = cd;
fileID= fopen(fullfile(BaseDir, 'HPC','Wally', 'WallyScript_Template.sh'));
[jobArrayTemplate] = textscan(fileID, '%s','Delimiter', '\n');
fclose(fileID);

%Replace job ID, Number of files and File path
JobArray = regexprep(jobArrayTemplate{1},'JOB_ID',JobID);

JobArray = regexprep(JobArray,'NUM_FILES',num2str(NumFrame));

JobArray = regexprep(JobArray,'VAR_FILE_PATH',VarListPath);

JobArray = regexprep(JobArray,'OUT_DIR',[VarListPath(1:end-11),'OUT']);



%Write Job array script
JobArrayPath = fullfile(NASfolderPath, ['batchArray_',JobID, '.sh']);
fileID= fopen(JobArrayPath, 'w+');

for L = 1:length(JobArray)
    fprintf(fileID,'%s\n',JobArray{L});
end
fclose(fileID);


