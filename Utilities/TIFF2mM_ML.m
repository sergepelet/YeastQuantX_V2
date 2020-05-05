%% Switching Data from Maoussi LHUILLIER-AKAKPO
%Main directory
ROOT = '/Volumes/DMF/GROUPS/gr_Pelet/Serge/DNAdamageSensor/Maoussi/161216';
%Folders to Move
FolderIN = {'161216_1432_2h30Gal', '161216_1432_4hGal'};
%{'161130_WT_4hGal', '161130_WTinc_4hGal','161130_mec1_4hGal','161130_mec1inc_4hGal','161130_sml1_4hGal','161130_sml1inc_4hGal'}; %,  '160412CDC5HA_23deg_4h',  '160412CDC5HA_32deg_1h' , '160412CDC5HA_32deg_2h' , '160412CDC5HA_32deg_4h';...
%   '160412CDC5HA_23deg_0h_', '160412CDC5HA_23deg_4h_', '160412CDC5HA_32deg_1h_', '160412CDC5HA_32deg_2h_', '160412CDC5HA_32deg_4h_'};

%Channel names
ChannelIN = {'BF', 'CFP', 'GFP', 'mCherry'};
Zstack = 15;

%Folder which will recieve the files
FolderOUT = '161216_GalTime';

ChannelOUT = {'01-BF', '02-CFP', '03-YFP', '04-RFP'};

%% Transfer for multiple position
mkdir(ROOT,FolderOUT)

for F = 1:length(FolderIN)
    PosFolderIN = fullfile(ROOT, FolderIN{F});
    FolderF = fullfile(ROOT, FolderIN{F},'*.TIF');
    
    ListTIFFiles= dir(FolderF);
    NumTIFfiles = length(ListTIFFiles);
    NumPosForCondition = round(length(ListTIFFiles)/length(ChannelIN));
    for P = 1:NumTIFfiles
        TFname = ListTIFFiles(P).name
        US = findstr(ListTIFFiles(P).name, '_');
        NumTIFPos = str2num(ListTIFFiles(P).name(US(3)+1:US(4)-1));
        
        PosFolderOUT = fullfile(ROOT,FolderOUT , ['Pos',num2str((F-1)*NumPosForCondition + NumTIFPos,'%02d')]);
        if ~exist(PosFolderOUT,'file')
            mkdir(PosFolderOUT)
        end
        
        %BF images save two Z frames
        if ~isempty(strfind(ListTIFFiles(P).name, ChannelIN{1}))
            ImgNameIN = fullfile(PosFolderIN, ListTIFFiles(P).name);
            ImgNameOUT = fullfile(PosFolderOUT, ['img_000000000_', ChannelOUT{1}, '0_000.tif']);
            Img = imread(ImgNameIN, 'TIFF', 4);
            imwrite(Img, ImgNameOUT, 'TIFF');
            ImgNameOUT = fullfile(PosFolderOUT, ['img_000000000_', ChannelOUT{1}, '1_000.tif']);
            Img = imread(ImgNameIN, 'TIFF', 10);
            imwrite(Img, ImgNameOUT, 'TIFF');
        end
        
        % Other Fluo Im
        for C = 2:length(ChannelIN)
            if ~isempty(strfind(ListTIFFiles(P).name, ChannelIN{C}))
                ImgNameIN = fullfile(PosFolderIN, ListTIFFiles(P).name);
                ImgNameOUT = fullfile(PosFolderOUT, ['img_000000000_', ChannelOUT{C}, '_000.tif']);
                for Z = 1:Zstack
                    Img = imread(ImgNameIN, 'TIFF', Z);
                    if Z==1
                        Zproject =Img;
                    else
                        Zproject = max(Zproject, Img);
                    end
                end
                imwrite(Zproject, ImgNameOUT, 'TIFF');
                
            end
        end
         Transfer{(F-1)*NumPosForCondition + NumTIFPos+1,1} = [PosFolderIN];
         Transfer{(F-1)*NumPosForCondition + NumTIFPos+1,2} = [PosFolderOUT];
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



