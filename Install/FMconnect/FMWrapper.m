function FMWrapper(RecordNumStr)



Commas = strfind(RecordNumStr, ',');
Commas = [1,Commas,length(RecordNumStr)];
for c = 1:length(Commas)-1
    SearchNum(c) = str2num(RecordNumStr(Commas(c)+1:Commas(c+1)-1));
end

button = questdlg('Select the platform for the analysis','Platform','Mac','PC', 'PC_parallel','Mac');

if strcmp(button, 'PC_parallel')
    System = 'PC_para';
else
    System = button;
end

cd ..
cd ..

VarCell = PrepareAnalysis('ExpNum', SearchNum, System);


button = questdlg('Run the analysis','Run','Yes','No','No');

if strcmp(button, 'Yes')
    ParallelImgAnalysis(VarCell)
end

quit