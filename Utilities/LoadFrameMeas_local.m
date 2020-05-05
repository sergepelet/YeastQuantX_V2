function  Var = LoadFrameMeas_local(Var)
tic

%Loop through all frames
for F = 1:Var.Analysis.NumFrame
    %Assemble Filepath to Frame.mat files
     OutPath = strrep(Var.Analysis.OutPath, '/groups/dmf/', '/Volumes/DMF/GROUPS/');
     %   OutPath = strrep(Var.Analysis.OutPath, '/mnt/spelet1/', '/Volumes/DMF/GROUPS/gr_Pelet/');
   Var.Analysis.OutPath = OutPath;
   ImgFolder = fileparts(Var.Analysis.OutPath);
    %ImgFolder = '/Users/serge/Projects/ImageAnalysis/MeasAnalysis/SungSik/Ouputfiles';
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(F,'%05d'),'.mat'])
    %Load file
    load(FileName)
    
    %For first frame create VarOUT
    if F == 1
        VarOUT = Var;
        VarOUT.SegMeas(F) = Var.Measurements;
        VarOUT = rmfield(VarOUT, 'Measurements');
    else  %Other frame add all measurements to SegMeas structure
        VarOUT.SegMeas(F) = Var.Measurements;
        VarOUT.Analysis.TimeStamp(F) = Var.Analysis.TimeStamp;
    end
end

%Set Var to VarOUT
Var = VarOUT;
        
%Save Timing Info
Var.Analysis.Timing.(mfilename) = toc;