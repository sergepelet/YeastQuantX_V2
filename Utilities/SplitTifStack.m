FolderList = {'/Volumes/DMF/GROUPS/gr_Pelet/Min_B/2017/170815/Hela_ERK_DS-ND_2', '/Volumes/DMF/GROUPS/gr_Pelet/Min_B/2017/170815/Hela_ERK_DS-ND_3','/Volumes/DMF/GROUPS/gr_Pelet/Min_B/2017/170815/Hela_ERK_DS-ND_4'};
%Illumination settings name in the order they are saved in the file
%IllumList = {'00-BF0','06-CFPtriple','07-YFPtriple','08-RFPtriple',  '00-BF1'};
IllumList = {'00-BF0','05-RFPdual','09-DAPIquad'};

%loop through all folders
for F = 1:length(FolderList)
    %Get list of file in folder
    FileList = dir(FolderList{F});
    %loop through all files
    for f = 1:length(FileList)
        %Get filename and extension
        [~,FileName, FileExt] = fileparts(FileList(f).name);
        %Verify is it is a tif file
        if FileList(f).isdir == 0 && strcmpi(FileExt, '.tif')
            %Find Pos number
            US = strfind(FileName, '_');
            PosStr = FileName(US(end)+1:end-4);
            %Make Pos folder
            PosFolder = fullfile(FolderList{F}, [PosStr]);
            mkdir(PosFolder);
            
            %Initalize Values
            Tpts = 0;
            Illum = 1;
            %Create Metadat file
            fileID = fopen(fullfile(PosFolder,'metadata.txt'), 'w+','n','UTF-8');
            
            %Read file info
            ImgFile = fullfile(FolderList{F}, FileList(f).name)
            ImgInfo = imfinfo(ImgFile);
            for N = 1:length(ImgInfo)
                %Read image from file
                Img = imread(ImgFile,'Index', N, 'Info', ImgInfo);
                %create single img file
                ImgName = ['img_', num2str(Tpts, '%09d'),'_',IllumList{Illum},'_000.tif'];
                imwrite(Img, fullfile(PosFolder,ImgName), 'tif');
                fprintf(fileID,'%s', ['FrameKey-',num2str(Tpts),'-', num2str(Illum-1),'-0": { "']);
                fprintf(fileID, '\n');
                fprintf(fileID,'%s', ['"FileName": "', ImgName, '"']);
                fprintf(fileID, '\n');
                MetaData = regexprep(ImgInfo(N).UnknownTags(end).Value, ',', ',\n');
                fprintf(fileID,'%s', MetaData);
                fprintf(fileID, '\n');
                %Increment Illum andTime point
                Illum = Illum+1;
                if Illum >length(IllumList)
                    Illum = 1;
                    Tpts = Tpts+1
                end
            end
            fclose(fileID);
        end
    end
end
