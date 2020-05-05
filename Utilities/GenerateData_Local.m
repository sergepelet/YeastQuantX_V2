
addpath('/Users/serge/Projects/ImageAnalysis/YeastQuantX/CommonProcess')

BFList  = {...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP928_0p2MNaCl_A_1';...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP928_0p2MNaCl_B_1';...
 '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP929_0p2MNaCl_A_1';...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP929_0p2MNaCl_B_1';...
     '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP930_0p2MNaCl_A_1';...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_ySP930_0p2MNaCl_B_1';...
     '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_yVW403_0p2MNaCl_A_1';...
    '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2020/200117/0117_yVW403_0p2MNaCl_B_1';...
    };
for F = 1:numel(BFList)

BaseFolder = BFList{F}

PosList = dir(fullfile(BaseFolder, 'Pos*'));
%% Loop through all MAT folders


for M = 1:length(PosList)
    
    DataFileName = fullfile(BaseFolder,'DATAOut', ['Pos_',num2str(M, '%03d'),'_Data.mat']);
    if exist(DataFileName,'file') == 2
        fprintf('\n DATA.mat exists\n')
        if M == 1
            load(DataFileName)
        end
    else
        DataOutFolder = fullfile(BaseFolder,'DATAOut');
         if ~exist(DataOutFolder,'file')
            fprintf('\n Create DATAOut\n')
            mkdir(DataOutFolder)
         end
        Frame1Path = fullfile(BaseFolder,['MATOut_', num2str(M, '%03d')], 'Frame_00001.mat');
        load(Frame1Path);
        %Load Frame.mat data
        Var = LoadFrameMeas_local(Var);
        assignin('base','Var_LoadFr',Var);
        %Track Main object
        Var.Analysis.TrackObj{1} = 'Nucl';
        
        Var = TrackObjectsSingle(Var, 1);
        
        assignin('base','VarTrack',Var);
        
        %Combine Frame measurements in single matrix
        Var = CombineFrame(Var);
        
        DataFileName = fullfile(BaseFolder,'DATAOut', ['Pos_',num2str(M, '%03d'),'_Data.mat'])
        % Save Var to Data.mat file
        save(DataFileName, 'Var', '-v7.3');
        close all
    end
end
Var = CombineData_local(Var, BaseFolder);

end
rmpath('/Users/serge/Projects/ImageAnalysis/YeastQuantX/CommonProcess')