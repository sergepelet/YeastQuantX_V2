function Var = FigureFlow(Var)

Var.Figure.Display = 'on';
%List of Subroutines which is used to define the figure number
iter = 0;
if isfield(Var.Analysis, 'LoadIllum')
    if iscell(Var.Analysis.LoadIllum)
        for L = 1:length(Var.Analysis.LoadIllum)
            iter = iter+1;
            Var.Figure.List{iter} = 'LoadImg' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'LoadImg' ;
    end
end

if isfield(Var.Analysis, 'Flatify')
    if iscell(Var.Analysis.LoadIllum)
        for L = 1:length(Var.Analysis.LoadIllum)
            iter = iter+1;
            Var.Figure.List{iter} = 'CorrFlat' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'CorrFlat' ;
    end
end

if isfield(Var.Analysis, 'SegPhiImg')
    if iscell(Var.Analysis.SegPhiImg)
        for L = 1:length(Var.Analysis.SegPhiImg)
            iter = iter+1;
            Var.Figure.List{iter} = 'SegPhi' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'SegPhi' ;
    end
end

if isfield(Var.Analysis, 'SegYFImg')
    if iscell(Var.Analysis.SegYFImg)
        for L = 1:length(Var.Analysis.SegYFImg)
            iter = iter+1;
            Var.Figure.List{iter} = 'SegYF' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'SegYF' ;
    end
    
end

if isfield(Var.Analysis, 'TrackObj')
    if iscell(Var.Analysis.TrackObj)
        for L = 1:length(Var.Analysis.TrackObj)
            iter = iter+1;
            Var.Figure.List{iter} = 'TrackObj' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'TrackObj' ;
    end
end


if isfield(Var.Analysis, 'SplitGroupIn')
    if iscell(Var.Analysis.SplitGroupIn)
        for L = 1:length(Var.Analysis.SplitGroupIn)
            iter = iter+1;
            Var.Figure.List{iter} = 'SplitCell' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'SplitCell' ;
    end
    
end
if isfield(Var.Analysis, 'AroundNucl')
    if iscell(Var.Analysis.AroundNucl)
        for L = 1:length(Var.Analysis.AroundNucl)
            iter = iter+1;
            Var.Figure.List{iter} = 'CellAroundNucl' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'CellAroundNucl' ;
    end
    
end


if isfield(Var.Analysis, 'ExpNuclNucl')
    if iscell(Var.Analysis.ExpNuclNucl)
        for L = 1:length(Var.Analysis.ExpNuclNucl)
            iter = iter+1;
            Var.Figure.List{iter} = 'ExpandNucl' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'ExpandNucl' ;
    end
    
end

if isfield(Var.Analysis, 'SecPrimObj')
    if iscell(Var.Analysis.SecPrimObj)
        for L = 1:length(Var.Analysis.SecPrimObj)
            iter = iter+1;
            Var.Figure.List{iter} = 'SecPrimObj' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'SecPrimObj' ;
    end
    
end
if isfield(Var.Analysis, 'ExpObjSmall')
    if iscell(Var.Analysis.ExpObjSmall)
        for L = 1:length(Var.Analysis.ExpObjSmall)
            iter = iter+1;
            Var.Figure.List{iter} = 'ExpObjSmall' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'ExpObjSmall' ;
    end
    
end


% 
% if isfield(Var.Analysis, 'CorrFlatIn')
%     if iscell(Var.Analysis.CorrFlatIn)
%         for L = 1:length(Var.Analysis.CorrFlatIn)
%             iter = iter+1;
%             Var.Figure.List{iter} = 'CorrFlat' ;
%         end
%     else
%         iter = iter+1;
%         Var.Figure.List{iter} = 'CorrFlat' ;
%     end
%     
% end


if isfield(Var.Analysis, 'SubtractIn')
    if iscell(Var.Analysis.SubtractIn)
        for L = 1:length(Var.Analysis.SubtractIn)
            iter = iter+1;
            Var.Figure.List{iter} = 'SubBackground' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'SubBackground' ;
    end
    
end


if isfield(Var.Analysis, 'RGBOut')
    if iscell(Var.Analysis.RGBOut)
        for L = 1:length(Var.Analysis.RGBOut)
            iter = iter+1;
            Var.Figure.List{iter} = 'RGBmix' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'RGBmix' ;
    end
    
end


if isfield(Var.Analysis, 'MovieImg')
    if iscell(Var.Analysis.MovieImg)
        for L = 1:length(Var.Analysis.MovieImg)
            iter = iter+1;
            Var.Figure.List{iter} = 'MakeMovie' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'MakeMovie' ;
    end
    
end

if isfield(Var.Analysis, 'DispMeasObject')
    if iscell(Var.Analysis.DispMeasObject)
        for L = 1:length(Var.Analysis.DispMeasObject)
            iter = iter+1;
            Var.Figure.List{iter} = 'DispMeas' ;
        end
    else
        iter = iter+1;
        Var.Figure.List{iter} = 'DispMeas' ;
    end
    
end
