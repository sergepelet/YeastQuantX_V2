BFList  = {'/Volumes/DMF/GROUPS/gr_Pelet/Serge/2019/190625/0625_Hog1YFP_Msn2RFP_A_1/';...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2019/190625/0625_Hog1YFP_Msn2RFP_B_1/';...
    };

for F = 1:numel(BFList)
    VarPath = fullfile(BFList{F}, 'DATAOut', 'Pos_001_Data.mat')
    load(VarPath)
    
    Var = CombineData_local(Var, BFList{F});
    
    AllFields = fieldnames(Export)
    
    for E = 1:length(AllFields)
        ExportALL.(AllFields{E}) = Export.(AllFields{E});
    end
end

