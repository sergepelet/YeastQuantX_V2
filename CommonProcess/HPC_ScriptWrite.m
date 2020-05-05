function HPC_ScriptWrite(JobID, NumFrame, VarListPath, NASfolderPath)
%Get VARLIST file info


%% Open template for jobarray script

BaseDir = cd;
fileID= fopen(fullfile(BaseDir, 'HPC', 'jobArray_VX_Template.sh'));
[jobArrayTemplate] = textscan(fileID, '%s','Delimiter', '\n');
fclose(fileID);

%Replace job ID, Number of files and File path
JobArray = regexprep(jobArrayTemplate{1},'JOB_ID',JobID);

JobArray = regexprep(JobArray,'NUM_FILES',num2str(NumFrame));

JobArray = regexprep(JobArray,'VAR_FILE_PATH',VarListPath);


%Write Job array script
JobArrayPath = fullfile(NASfolderPath, ['jobArray_',JobID, '.sh']);
fileID= fopen(JobArrayPath, 'w+');

for L = 1:length(JobArray)
    fprintf(fileID,'%s\n',JobArray{L});
end
fclose(fileID);
