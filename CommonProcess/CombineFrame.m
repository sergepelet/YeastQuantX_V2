function Var = CombineFrame(Var)
tic
%Get number cells in data set
NumberCell = max(cellfun(@max, Var.TrackLabel));
Var.Analysis.NumberCell = NumberCell;

%Get the Position number 
CurrentPos = Var.Analysis.CurrentPos;

%Get names of all objects and loop
ObjFields = fieldnames(Var.SegMeas(1));

for Obj = 1:length(ObjFields)
  % Obj= Obj
    %Get names of all illuminations 
    IllumNames = Var.Analysis.MeasureIntensity;
    %Check if illum were measured in object
    if isfield(Var.SegMeas(1).(ObjFields{Obj}), IllumNames{1})
        %If measured then loop
        for I = 1:length(IllumNames)
         %   I = I
            %Get the name of all measurements
           MeasFields = fieldnames(Var.SegMeas(1).(ObjFields{Obj}).(IllumNames{I}));
            %loop through all measurements
            for M = 1:length(MeasFields)
           %   M = M
                %Initalize measurement matrix filled with zeros
                if isnumeric(Var.SegMeas(1).(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}))
                    Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}) = zeros(Var.Analysis.NumFrame, NumberCell);
                end
%                 assignin('base', 'VarInit', Var )
%                 SM = size(Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}))
                % check if meas is checkLabel
                if strcmp(MeasFields{M}, 'CheckLabel')
                    %loop through all cells
                    for C =  1:NumberCell
                        %Calculate CheckLabel by assembling CellNumber and current position value divided by 1000
                        CheckLabel = C+CurrentPos/1000;
                        %Loop though all frames
                        for F = 1:Var.Analysis.NumFrame
                            %Check if cell is tracked at current frame
                            CellInd = find(Var.TrackLabel{F} == C);
                            if ~isempty(CellInd)
                            %    IndMat(F,C) = CellInd; 
                                %Set value of checklabel
                                Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M})(F,C) = CheckLabel;
                                
                                Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).VerifyLabel(F,C) = ...
                                    Var.SegMeas(F).(ObjFields{Obj}).SegLabel(CellInd);
                                
                                Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).CheckCheckLabel(F,C) = ...
                                    Var.SegMeas(F).(ObjFields{Obj}).(IllumNames{I}).CheckLabel(CellInd);
                                
                                %Set PixelList value if Illum  is one
                                if I == 1
                                    Var.Measurements.(ObjFields{Obj}).PixelList{F,C} = Var.SegMeas(F).(ObjFields{Obj}).PixelList{CellInd};
                                    Var.Analysis.TimeAxis = Var.Analysis.TimeStamp-Var.Analysis.TimeZero;
                                end
                            end
                        end
                    end
                else %For other measurements
                    if isnumeric(Var.SegMeas(1).(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}))
                        %loop through all cells
                        for C = 1:NumberCell
                           %   C = C
                            %Loop though all frames
                            for F = 1:Var.Analysis.NumFrame
                             %       F =F
                                %Get cell index from trackMatrix
                                CellInd = find(Var.TrackLabel{F} == C);
                                %Set value of measurment in Measurement matrix
                                if ~isempty(CellInd)
                                    Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M})(F,C) = ...
                                        Var.SegMeas(F).(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M})(CellInd);
                                end
                            end
                        end
                        
                    else
                            %loop through all cells
                        for C = 1:NumberCell
                            %       C = C
                            %Loop though all frames
                            for F = 1:Var.Analysis.NumFrame
                                %         F =F
                                %Get cell index from trackMatrix
                                CellInd = find(Var.TrackLabel{F} == C);
                                %Set value of measurment in Measurement matrix
                                if ~isempty(CellInd)
                                    Var.Measurements.(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}){F,C} = ...
                                        Var.SegMeas(F).(ObjFields{Obj}).(IllumNames{I}).(MeasFields{M}){CellInd};
                                end
                            end
                        end 
                    end
                end
            end
        end
    end
end

%assignin('base','IndMat',IndMat);
%Remove SegMeas structure
%Var = rmfield(Var, 'SegMeas');
%Save Timing Info
Var.Analysis.Timing.(mfilename) = toc;




