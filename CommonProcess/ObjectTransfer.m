function Var = ObjectTransfer(Var)
tic
%Transfer objects from previous segmentation to new frame when images for
%segmentation were skipped

%Assemble Filepath to Frame.mat files
    ImgFolder = fileparts(Var.Analysis.OutPath);
    FileName = fullfile(ImgFolder,['MATOut_', num2str(Var.Analysis.CurrentPos, '%03d')],['Frame_', num2str(Var.Analysis.SegmentedFrame,'%05d'),'.mat']);
    %Load file
    Prev = load(FileName)
 assignin('base','Prev',Prev);
%Get Object list
ObjectList = fieldnames(Var.Measurements);
%FI = Var.Analysis.FrameIter

%For all objects transfer centroid location, Pixel list and label
for O = 1:length(ObjectList)
   try
%    ObjectList{O} 
    Var.Measurements.(ObjectList{O}).CenterX = Prev.Var.Measurements.(ObjectList{O}).CenterX;
    Var.Measurements.(ObjectList{O}).CenterY = Prev.Var.Measurements.(ObjectList{O}).CenterY;
    Var.Measurements.(ObjectList{O}).SegLabel = Prev.Var.Measurements.(ObjectList{O}).SegLabel;
    Var.Measurements.(ObjectList{O}).PixelList = Prev.Var.Measurements.(ObjectList{O}).PixelList;
    end
end
%Save Timing Info
Var.Analysis.Timing.(mfilename) = toc;