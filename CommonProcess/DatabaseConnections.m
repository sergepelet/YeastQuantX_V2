function Var = DatabaseConnections(Var)
% %% settings
% Datasource = 'local_11';  % local_10 or FMserver local_11
% DatabaseName= 'YeastQuant';
% User = 'matlab';
% Passwd = 'yq_matlab';
% ServerName = 'localhost';  %'localhost:2399'


%%


% Set preferences For database with setdbprefs.
s.DataReturnFormat = 'cellarray';
s.ErrorHandling = 'store';
s.NullNumberRead = 'NaN';
s.NullNumberWrite = 'NaN';
s.NullStringRead = 'null';
s.NullStringWrite = 'null';
%s.UseRegistryForSources = 'yes';
%s.TempDirForRegistryOutput = '';
setdbprefs(s)

%%


if strcmp(Var.Database.Version, 'FM10') &&  strcmp(Var.Database.ServerPort, 'localhost')
    %Set Database specific variables:

    %Filemaker server: here local
    Var.Database.ServerPort = [ServerName, ':2399'];


    
    %%Make connection to database.  Using JDBC driver.
    connection = database(Var.Database.Name,Var.Database.User,Var.Database.Psswd,'com.ddtek.jdbc.sequelink.SequeLinkDriver',...
        ['jdbc:sequelink://',Var.Database.ServerPort,';serverDataSource=',Var.Database.Name,';user=', Var.Database.User,';password=',Var.Database.Psswd]);

elseif strcmp(Var.Database.Version, 'FM11')
   % FMstring = ['jdbc:filemaker://',Var.Database.ServerPort,'/',Var.Database.Name]
    %%Make connection to database.  Using JDBC driver.
    connection = database(Var.Database.Name,Var.Database.User,Var.Database.Psswd,'com.filemaker.jdbc.Driver',...
        ['jdbc:filemaker://',Var.Database.ServerPort,'/',Var.Database.Name]);

end
    
%Check Connections use different code based on version number
if verLessThan('matlab','9.3')
   TestConnection = isconnection(connection);
else
    TestConnection = isopen(connection);
end

if TestConnection
    Var.Database.connection = connection;
    % If the connection failed, print the error message
else
    display(sprintf('Connection failed: %s', connection.Message));

    error('Error while connecting to the database') 
end

%%Get Experimental Record from the databse
%Build SQL search query
%Build Select string for numbers
if isnumeric(Var.Database.SearchStr)
    SelectStr = ['SELECT * FROM ',Var.Database.Table ,' WHERE ', Var.Database.SearchField, ' = '];
    for L = 1:length(Var.Database.SearchStr)
        if L < length(Var.Database.SearchStr)
            SelectStr = [SelectStr, num2str(Var.Database.SearchStr(L)), ' OR ', Var.Database.SearchField, ' = '];
        else
            SelectStr = [SelectStr, num2str(Var.Database.SearchStr(L))];
        end
    end
else
%Build Select string for String or Cell array
    SelectStr = ['SELECT * FROM ',Var.Database.Table ,' WHERE ', Var.Database.SearchField, ' IN ('''];
    if iscell(Var.Database.SearchStr)
        for i = 1:length(Var.Database.SearchStr)
            if i == 1
                SelectStr = [SelectStr, Var.Database.SearchStr{i}];
            else
                SelectStr = [SelectStr, ''', ''', Var.Database.SearchStr{i}];
            end
        end
    else
        SelectStr = [SelectStr, Var.Database.SearchStr];
    end
    SelectStr = [SelectStr, ''')'];
end
        
%Perform database search
Cursor = exec(Var.Database.connection, SelectStr);
Cursor = fetch(Cursor);
%Get data
Var.Database.Data = Cursor.Data;
%Get fieldnames
%With newer Matlab versions:
% Var.Database.FieldString = columnnames(Cursor,1);
%Safer solution:
%Get string with all names
SingleString = columnnames(Cursor);
%Find limits of column names
Delimiter = strfind(SingleString,',');
Delimiter = [0, Delimiter, length(SingleString)+1];
%Remove FieldString if existing
try
    Var.Database = rmfield(Var.Database, 'FieldString');
end
%Split Sting in a cell array
for i = 1:length(Delimiter)-1
    Var.Database.FieldString{i} = SingleString(Delimiter(i)+2:Delimiter(i+1)-2);
end
%Get number of records found
Var.Database.NumFoundRecord = size(Cursor.Data, 1);

if Var.Database.NumFoundRecord == 0
    error('No record found for:', SelectStr)
end


close(Cursor)
% Close database connection.
close(Var.Database.connection)

Var.Database = rmfield(Var.Database, 'connection');
