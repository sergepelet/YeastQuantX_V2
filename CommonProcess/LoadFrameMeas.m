function  Var = LoadFrameMeas(Var)
tic

%Loop through all frames
for F = 1:Var.Analysis.NumFrame
    %Assemble Filepath to Frame.mat files
    ImgFolder = fileparts(Var.Analysis.OutPath);
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(F,'%05d'),'.mat']);
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
    VarOUT.Analysis.FrameTiming(F) = Var.Analysis.Timing;
end

%Set Var to VarOUT
Var = VarOUT;
        
%Save Timing Info
Var.Analysis.Timing.(mfilename) = toc;