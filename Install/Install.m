%% Install Script for YeastQUant
%Works for Filemaker 11
%Copies fmjdbc.jar to the Install/FMconnect folder
%Add fmjdbc.jar to the java classpath of matlab
%Set-up the connection with the filemaker database
%Test the connection
%Enter the image server name
%On the Mac edit a script to have direct link from filemaker to matlab

%% Get folder names
InstallFolder = cd;
cd ..
MatlabCodeFolder = cd;
cd(InstallFolder)
FMconnectFolder = fullfile(InstallFolder, 'FMconnect');

%% Locate fmjdbc.jar driver file in install CD from Filemaker and copy to FMconnect folder

%check if drivers are already in FMconnect folder

if exist(fullfile(FMconnectFolder, 'fmjdbc.jar'), 'file') ~= 2
    
    [FmDriver,FmPathName] = uigetfile('*.jar','Locate fmjdbc.jar file on FileMake install disc');
    
    if ~strcmp(FmDriver, 'fmjdbc.jar')
        
        helpdlg('The fmjdbc.jar file can be found on the Install disc of filemaker at \FileMaker Pro 11 Advanced\xDBC\JDBC Client Driver Installer','fmjdbc File Help')
        pause(2)
        [FmDriver,FmPathName] = uigetfile('*.jar','Locate fmjdbc.jar file on FileMake install disc');
        
        if ~strcmp(FmDriver, 'fmjdbc.jar')
            msgbox('Instal process stopped: Selected .jar file is not the correct one','filemaker driver selection','error')
            
            return
        end
    end
    %Copy Jar file to folder
    copyfile(fullfile(FmPathName,FmDriver),fullfile(FMconnectFolder,FmDriver))

else
    button = questdlg('Reload JDBC drivers for Filemaker to matlab connection ','Filemaker Driver', 'No');
    if strcmp(button, 'Yes')
        [FmDriver,FmPathName] = uigetfile('*.jar','Locate fmjdbc.jar file on FileMake install disc');
        %Copy Jar file to folder
        copyfile(fullfile(FmPathName,FmDriver),fullfile(FMconnectFolder,FmDriver))

    elseif strcmp(button, 'No')
        
        FmDriver = 'fmjdbc.jar';
        
    elseif strcmp(button, 'Cancel')
        return
    end
    
end



%% Add FMconnect to dynamic path of Matlab


javaaddpath(fullfile(FMconnectFolder,FmDriver))


%% Set-up of Matlab to database connection

%Get Filemaker version
% FMVersion = listdlg('PromptString','Choose FileMaker version',...
%                 'SelectionMode','single',...
%                 'ListString',{'FileMaker 10','FileMaker 11'}, 'ListSize',[150,75]);
% if FMVersion == 1
%     Database.Version = 'FM10';
% elseif FMVersion == 2
%    Database.Version = 'FM11';
%end


if exist( fullfile(MatlabCodeFolder,'CommonProcess','DatabaseConnextionParameter.mat'), 'file') == 2
    load(fullfile(MatlabCodeFolder,'CommonProcess','DatabaseConnextionParameter.mat'))
    DlgText = {'Enter database name','Enter username','Enter password', 'Enter server adress'};
    DlgTitle = 'Database connection parameters';
    DlgDefault = {Database.Name,Database.User,Database.Psswd, Database.ServerPort};
    if strcmpi(Database.ServerPort, 'localhost')
        Platform = 1;
    else
        Platform = 2;
    end
else
    
    Database.Version = 'FM11';
    
    
    %Get Server or local
    Platform = listdlg('PromptString','Choose database type',...
        'SelectionMode','single',...
        'ListString',{'Local','Server'}, 'ListSize',[150,75]);
    
    %Get database parameters
    if Platform == 1
        DlgText = {'Enter database name','Enter username','Enter password'};
        DlgTitle = 'Database connection parameters';
        DlgDefault = {'YeastQuant','matlab','yq_matlab'};
    elseif Platform == 2
        DlgText = {'Enter database name','Enter username','Enter password', 'Enter server adress'};
        DlgTitle = 'Database connection parameters';
        DlgDefault = {'YeastQuant','matlab','yq_matlab', 'filemaker.domain.ch'};
    end
end


DlgAnswer = inputdlg(DlgText,DlgTitle,1,DlgDefault);

Database.Name = DlgAnswer{1};
Database.User = DlgAnswer{2};
Database.Psswd = DlgAnswer{3};
if length(DlgAnswer) == 3
    Database.ServerPort = 'localhost'; 
else
    Database.ServerPort = DlgAnswer{4};
end

%For Filemaker 10
%
%      %%Make connection to database.  Using JDBC driver.
%     connection = database(Database.Name,Database.User,Database.Psswd,'com.ddtek.jdbc.sequelink.SequeLinkDriver',...
%         ['jdbc:sequelink://',Database.ServerPort,';serverDataSource=',Database.Name,';user=', Database.User,';password=',Database.Psswd]);

%For Filemaker 11 and 12

%%Make connection to database.  Using JDBC driver.
connection = database(Database.Name,Database.User,Database.Psswd,'com.filemaker.jdbc.Driver',...
    ['jdbc:filemaker://',Database.ServerPort,'/',Database.Name]);


%Test Connections
if isopen(connection)
    
    msgbox(['Successful connection to ', Database.Name],'Database Connection', 'modal' )
    pause(2)
    % Close database connection.
    close(connection)
    
    %Add fm driver path to Database structure
	Database.DriverPath = fullfile(FMconnectFolder,FmDriver);
    %Save database settings in install folder
    save(fullfile(MatlabCodeFolder,'CommonProcess','DatabaseConnextionParameter.mat'), 'Database')
    
    % If the connection failed, print the error message
else
    if Platform == 1
        msgbox(['Could not connect to ', Database.Name, ' Verify that the database is opened locally and that the connection parameters are correct: ',connection.Message],'Database Connection','error', 'modal')
    elseif Platform == 2
        msgbox(['Could not connect to ', Database.Name, ' Verify that the database is properly shared on your server and that the connection parameters are correct: ',connection.Message],'Database Connection','error','modal')
    end
    display(sprintf('Connection failed: %s', connection.Message));
    return
end


%% Specify image server Name

if exist( fullfile(MatlabCodeFolder,'CommonProcess','ImageStoragePath.mat'), 'file') == 2
    load(fullfile(MatlabCodeFolder,'CommonProcess','ImageStoragePath.mat'))
    DlgText = {'Enter Mac server path','Enter PC server path'};
    DlgTitle = 'Server Path';
    DlgDefault = {ImageStoragePath.MAC,ImageStoragePath.PC};
    DlgOptions.Resize = 'on';
    
    DlgAnswer = inputdlg(DlgText,DlgTitle,1,DlgDefault, DlgOptions);
    ImageStoragePath.MAC = DlgAnswer{1};
    ImageStoragePath.PC = DlgAnswer{2};
else
    
    ImgServer = listdlg('PromptString','Choose Image source',...
        'SelectionMode','single',...
        'ListString',{'Images stored locally','Images on a server'}, 'ListSize',[150,75]);
    if ImgServer == 1
        ImageStoragePath.MAC = 'local';
        ImageStoragePath.PC = 'local';
    elseif ImgServer == 2
        DlgText = {'Enter Mac server path','Enter PC server path'};
        DlgTitle = 'Server Path';
        DlgDefault = {'/Volumes/ImgServerName/','\\nas.domain.ch\ImgServerName\'};
        DlgOptions.Resize = 'on';
        
        DlgAnswer = inputdlg(DlgText,DlgTitle,1,DlgDefault, DlgOptions);
        ImageStoragePath.MAC = DlgAnswer{1};
        ImageStoragePath.PC = DlgAnswer{2};
        
    end
end
%Save Server settings to CommonProcesses folder

save(fullfile(MatlabCodeFolder,'CommonProcess','ImageStoragePath.mat'), 'ImageStoragePath')


%% Add direct activation of Matlab from FileMaker

% works only on MAC because it requires Applescript from filemaker to start
% matlab

if ismac
    
    
    %% Edit YQlink.txt File
    
    %Open the YQlink.txt file
    fileID = fopen(fullfile(FMconnectFolder,'YQlink.txt'), 'rt+','n','UTF-8');
    frewind(fileID);
    YQtxt = textscan(fileID, '%s', 'Delimiter', '\t');
    fclose(fileID);
    
    
    %Replace _YQCodeFolder_ with actual location of Matlab code
    Limit = strfind(YQtxt{1}{1}, '_');
    YQCommand{1} = [YQtxt{1}{1}(1:Limit(1)-1),MatlabCodeFolder, YQtxt{1}{1}(Limit(2)+1:end)];
    
    %Replace MATLABROOT with actual location of Matlab instalation
    Limit = strfind(YQtxt{1}{2}, '_');
    YQCommand{2} =[YQtxt{1}{2}(1:Limit(1)-1),matlabroot, YQtxt{1}{2}(Limit(2)+1:end)];
    
    
    %create YQlink.command file
    YQName = ['Link_', Database.Name ,'.command'];
    fileID = fopen(fullfile(FMconnectFolder,YQName), 'wt+','n','UTF-8');
    %Write to file
    fprintf(fileID,'%s', YQCommand{1});
    fprintf(fileID, '\n');
    fprintf(fileID,'%s', YQCommand{2});
    fclose(fileID);
    
    %Set file as executable
    fileattrib(fullfile(FMconnectFolder,YQName),'+x')
    %Creat FMconnect folder in Library of user and copy YQlimk.command to
    %this folder
    mkdir('~/Library', 'FMconnect')
    movefile(fullfile(FMconnectFolder,YQName),fullfile('~/Library/FMconnect/',YQName))
    
    
end

%% Finish instalation message

msgbox('Successful Installation','Install Process', 'modal' )
pause(3)
%clear all variables
clear all

