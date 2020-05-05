function Var = RemoveMeasurements(Var)

ObjectList = Var.Analysis.MeasureObject;
IllumList = Var.Analysis.MeasureIntensity;

for O = 1:length(ObjectList)
    for I = 1:length(IllumList)
       Var.Measurements.(ObjectList{O}) = rmfield(Var.Measurements.(ObjectList{O}), IllumList{I});
        
    end
end
