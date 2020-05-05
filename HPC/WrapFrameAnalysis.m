function WrapFrameAnalysis(FileName)

warning off

load(FileName)
fprintf('Loaded Var: ');
VarFrame.Figure.Display = 'off';

ImgFolder = fileparts(VarFrame.Analysis.OutPath);
	fprintf([ImgFolder, ' Pos: ', num2str(VarFrame.Analysis.CurrentPos), ' Frame: ' num2str(VarFrame.Analysis.CurrentFrame),'\n'])
	MATOutFolder = fullfile(ImgFolder,['MATOut_', num2str(VarFrame.Analysis.CurrentPos, '%03d')]);
if exist(MATOutFolder, 'file') == 0
	mkdir(MATOutFolder)
end


DataOutFolder = fullfile(ImgFolder,'DATAOut');
if exist(DataOutFolder, 'file') == 0
     mkdir(DataOutFolder)
end

FileName = fullfile(ImgFolder,['MATOut_', num2str(VarFrame.Analysis.CurrentPos, '%03d')],['Frame_', num2str(VarFrame.Analysis.CurrentFrame,'%05d'),'.mat']);

if exist(FileName, 'file') == 0
	FrameAnalysis(VarFrame);
else
	fprintf('MAT file already analyzed\n')
end
warning on
quit