function Var = GetDBFields(Var, CallNum)
%Get fieldnames and field content based on Var.Database obtained from Database Connection
%Get variable from structure
% Raw cell array from SQL output
DBData = Var.Database.Data;
%Fieldnames from SQL
FieldString = Var.Database.FieldString;
%Table name used to define variable structure
Table = Var.Database.Table;

%Check if Table is already a field in Var
if isfield(Var, Table)
    %If yes add new entry in structure
    NumTable = length(Var.(Table))+1;
else
    NumTable = 1;
end

%Distribute Database entries in Var
for F = 1:length(FieldString)
   % F = F
    %Check if entries are not null or Nan
    if ~(strcmp(DBData{CallNum,F}, 'null') | isnan(DBData{CallNum,F}))
        %Check if fieldstring contain numeric value
        if sum(isstrprop(FieldString{F}, 'digit')) == 0
            Var.(Table)(NumTable).(FieldString{F}) = DBData{CallNum,F};
        else
            %If fieldstring contains numeric value
            %Need to distribute data in same fieldname
            %Check where underscore is and how many there are
            us = strfind(FieldString{F}, '_');
            if length(us) == 1
            	%Get number of fields present in the database
                NCS = [FieldString{F}(us+1:end)];
                NumCol = str2num(NCS);
                Var.(Table)(NumTable).(FieldString{F}(1:us-1)){NumCol} = DBData{CallNum,F};
            elseif length(us) == 2
            %If 2 underscore => cell matrix
            %Define number of columns and rows
                NCS = [FieldString{F}(us(2)+1:end)];
                NumCol = str2num(NCS);
                NLS = [FieldString{F}(us(1)+1:us(2)-1)];
                NumLine = str2num(NLS);
                Var.(Table)(NumTable).(FieldString{F}(1:us-1)){NumLine, NumCol} = DBData{CallNum,F};
            end
        end
    end
end