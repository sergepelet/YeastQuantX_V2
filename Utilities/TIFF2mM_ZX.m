%% Switching Data from Zhou Xu
%Main directory
ROOT = '/Volumes/DMF/GROUPS/gr_Pelet/Serge/DNAdamageSensor/cdc13-1_zx';
%Folders to Move
FolderIN = {'160412CDC5HA_23deg_0h',  '160412CDC5HA_23deg_4h',  '160412CDC5HA_32deg_1h' , '160412CDC5HA_32deg_2h' , '160412CDC5HA_32deg_4h';...
            '160412CDC5HA_23deg_0h_', '160412CDC5HA_23deg_4h_', '160412CDC5HA_32deg_1h_', '160412CDC5HA_32deg_2h_', '160412CDC5HA_32deg_4h_'};

%Channel names
ChannelIN = {'1_ORG', '2_ORG', '3_ORG'};

%Folder which will recieve the files
FolderOUT = '160412_ZX_yT976_CDC13ts';

ChannelOUT = {'01-Phase', '02-CFP', '03-RFP'};

%% Transfer for multiple position
mkdir(ROOT,FolderOUT)
NumPos = 0;
for F = 1:length(FolderIN)
    %FolderF = fullfile(ROOT, [FolderIN{F},filesep,'Pos*']);
    if size(FolderIN,1) == 1
        FolderF = fullfile(ROOT, [FolderIN{F},'*']);
    else
        FolderF = fullfile(ROOT, FolderIN{1,F}, [FolderIN{2,F},'*']);
    end
       
    ListPosFolders= dir(FolderF);
    
    for P = 1:length(ListPosFolders)
        if isstrprop(ListPosFolders(P).name(length(FolderIN{end,F})+1), 'digit')
            NumPos = NumPos +1;
            if size(FolderIN,1) == 1
                PosFolderIN = fullfile(ROOT, ListPosFolders(P).name)
            else
                PosFolderIN = fullfile(ROOT, FolderIN{1,F}, ListPosFolders(P).name)
            end
            PosFolderOUT = fullfile(ROOT,FolderOUT , ['Pos',num2str(NumPos,'%02d')])
            mkdir(PosFolderOUT)
            
            ListImage= dir(PosFolderIN);
            for Img = 1:length(ListImage)
                for C = 1:length(ChannelIN)
                    if ~isempty(strfind(ListImage(Img).name, ChannelIN{C}))
                        ImgNameIN = fullfile(PosFolderIN, ListImage(Img).name);
                        ImgNameOUT = fullfile(PosFolderOUT, ['img_000000000_', ChannelOUT{C}, '_000.tif']);
                         movefile(ImgNameIN,ImgNameOUT)
                    end
                end
            end
        end

        
              
        %
                 Transfer{NumPos,1} = [PosFolderIN];
                 Transfer{NumPos,2} = [PosFolderOUT];
    end
end
save(fullfile(ROOT,FolderOUT,'Transfer.mat'), 'Transfer')


%% OTHER TRANSFER with single position
% 
% 
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

        
    
    