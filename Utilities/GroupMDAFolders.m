clear all
%Main directory
ROOT = '/Volumes/DMF/GROUPS/gr_Pelet/Serge/2019/190327';
%Folders to Move
FolderIN = {'0327_A_cyc_1','0327_A_arr_1','0327_A_15_1','0327_A_30_1','0327_A_45_1','0327_A_60_1','0327_A_75_1','0327_A_90_1','0327_A_105_1', '0327_A_120_1'}; %, '0621_SqV_Kss1nls_log_1', '0621_SqV_Kss1nls_aF_1',...
    %'0621_SqV_kss1d_log_1', '0621_SqV_kss1d_aF_1',...
   %  '0621_SqV_f3d_Kss1wt_log_1','0621_SqV_f3d_Kss1wt_aF_1', '0621_SqV_f3d_Kss1nls_log_1','0621_SqV_f3d_Kss1nls_aF_1'};


%Folder which will recieve the files
FolderOUT = '0327_LineA';

%% Transfer for multiple position in each folder
mkdir(ROOT,FolderOUT)
NumPos = 0;
for F = 1:length(FolderIN)
    %FolderF = fullfile(ROOT, [FolderIN{F},filesep,'Pos*']);
     FolderF = fullfile(ROOT, [FolderIN{F},filesep,'0*']);
    ListPosFolders= dir(FolderF);
    for P = 1:length(ListPosFolders)
        NumPos = NumPos +1;
        
        PosFolderIN = fullfile(ROOT, FolderIN{F}, ListPosFolders(P).name)
        PosFolderOUT = fullfile(ROOT,FolderOUT , ['Pos',num2str(NumPos,'%02d')])
        

       movefile(PosFolderIN,PosFolderOUT)
        
        Transfer{NumPos,1} = [FolderIN{F}, filesep, ListPosFolders(P).name];
         Transfer{NumPos,2} = [FolderOUT, filesep,'Pos',num2str(NumPos,'%02d')];
    end
end
save(fullfile(ROOT,FolderOUT,'Transfer.mat'), 'Transfer')


%% OTHER TRANSFER with single position


% mkdir(ROOT,FolderOUT)
% NumPos = 0;
% for F = 1:length(FolderIN)
%     
%         NumPos = NumPos +1;
%         
%         PosFolderIN = fullfile(ROOT, FolderIN{F});
%         PosFolderOUT = fullfile(ROOT,FolderOUT , ['Pos',num2str(NumPos,'%02d')]);
%         
%         movefile(PosFolderIN,PosFolderOUT)
%         
%         Transfer{NumPos,1} = [FolderIN{F}];
%          Transfer{NumPos,2} = [FolderOUT, filesep,'Pos',num2str(NumPos,'%02d')];
% end
% save(fullfile(ROOT,FolderOUT,'Transfer.mat'), 'Transfer')
% 
%         
%     FolderList = dir(ROOT);
%     mkdir(ROOT,FolderOUT)
% NumPos = 0;
% for F = 3:length(FolderList)
%     
%         NumPos = NumPos +1;
%         
%         %PosFolderIN = fullfile(ROOT, FolderIN{F});
%         PosFolderIN = fullfile(ROOT, FolderList(F).name);
%         PosFolderOUT = fullfile(ROOT,FolderOUT , ['Pos',num2str(NumPos,'%02d')]);
%         
%         movefile(PosFolderIN,PosFolderOUT)
%         
%         Transfer{NumPos,1} = [FolderList(F).name];
%          Transfer{NumPos,2} = [FolderOUT, filesep,'Pos',num2str(NumPos,'%02d')];
% end
% save(fullfile(ROOT,FolderOUT,'Transfer.mat'), 'Transfer')

        
