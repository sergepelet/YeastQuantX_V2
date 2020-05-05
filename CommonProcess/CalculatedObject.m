function Var = CalculatedObject(Var, CallNum)
tic

if nargin == 1
    CallNum = 1;
end

Debug =0;

%Assign Objects and intensities name to variables


FirstInt = Var.Analysis.CalcFirstInt{CallNum};
SecInt = Var.Analysis.CalcSecInt{CallNum};
FirstObj = Var.Analysis.CalcFirstObj{CallNum};
SecObj = Var.Analysis.CalcSecObj{CallNum};

Operation = Var.Analysis.CalcOperation{CallNum};

ObjOUT = Var.Analysis.CalcObjOUT{CallNum};

CurrentFrame = Var.Analysis.CurrentFrame;
%Transfer Geometrical properties from object to new Object

% Var.Measurements.(ObjOUT).CenterX{FrameIter} = Var.Measurements.(FirstObj).CenterX{FrameIter};
% Var.Measurements.(ObjOUT).CenterY{FrameIter} = Var.Measurements.(FirstObj).CenterY{FrameIter};
Var.Measurements.(ObjOUT).PixelList = Var.Measurements.(FirstObj).PixelList;
Var.Measurements.(ObjOUT).SegLabel = Var.Measurements.(FirstObj).SegLabel;

%Get fieldnames from measurements
MeasNames = fieldnames(Var.Measurements.(FirstObj).(FirstInt));
ExcludeList = {'CenterX','CenterY','CheckLabel','Area', 'Diameter', 'HiPixList', 'LoPixList', 'ConnectedHiPixList'}; %,'Polarity'};

%Loop through all measurements
%To make the new calculation
for M = 1:length(MeasNames)
    %If measurement is in exculde list Take the same as in first object
    %first int
    if max(strcmp(MeasNames{M}, ExcludeList)) == 1
        Var.Measurements.(ObjOUT).(FirstInt).(MeasNames{M}) = Var.Measurements.(FirstObj).(FirstInt).(MeasNames{M});
    else
        %Calculate difference between first meas and second measurement
        if strcmp(Operation, 'Difference')
            
%             Var.Measurements.(ObjOUT).(FirstInt).(MeasNames{M}) = cellfun(@(x,y) x-y,... 
%                 Var.Measurements.(FirstObj).(FirstInt).(MeasNames{M}),...
%                 Var.Measurements.(SecObj).(SecInt).(MeasNames{M}), 'UniformOutput', false);
%             
              Var.Measurements.(ObjOUT).(FirstInt).(MeasNames{M}) = Var.Measurements.(FirstObj).(FirstInt).(MeasNames{M}) - ...
                Var.Measurements.(SecObj).(SecInt).(MeasNames{M});
            
        elseif strcmp(Operation, 'Ratio')
             Var.Measurements.(ObjOUT).(FirstInt).(MeasNames{M}) = Var.Measurements.(FirstObj).(FirstInt).(MeasNames{M}) ./ ...
                Var.Measurements.(SecObj).(SecInt).(MeasNames{M});
            
        end
    end
end

    
%Save Timing Info
Var.Analysis.Timing.(mfilename)(CallNum) = toc;