function tablePrinter(header, csvData, tableTitle, folderName)

% Prints to csv file with header on first row.

%

%function call: tablePrinter(header, csvData, tableTitle)

%    header and tableTitle should be strings, csvData is numeric.



%go to the output directory

currentDirectory = pwd;

cd(folderName);



%make the output table

tableName = strcat(tableTitle,'.csv');

outid = fopen(tableName, 'w+');

fprintf(outid, '%s', header);

dlmwrite (tableName,csvData,'roffset',1,'-append','precision', 6);

fclose(outid);


%get back out to the original directory

cd(currentDirectory);



