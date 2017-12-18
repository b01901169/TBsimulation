function stataOutputMaker(recNo, folderName)
%building a stata output generator for the TB simulation

%grab recNo 
currentDirectory = pwd;
cd(folderName);
load('recNo.mat');
cd(currentDirectory);

%get rid of all-dead or all-not-yet-born people
allDead = false(size(intRecNo,1), 1);
healthMat = floor( (  mod(intRecNo,48) )/6 );
for rowNum = 1:size(intRecNo, 1)
    %if healthMat for that person is all 5 or all 6
    allDead(rowNum, 1) = (any(healthMat(rowNum,:) ~= 5) | any(healthMat(rowNum,:) ~= 6) );
end

%grab only full lifetime people
fullLifeMat = false(size(intRecNo,1), 1);
for rowNum = 1:size(intRecNo, 1)
    %if healthMat for that person has 5 and has 6
    fullLifeMat(rowNum, 1) = (any(healthMat(rowNum,:) == 6) && any(healthMat(rowNum,:) == 5) );
end

finalKeepers = allDead & fullLifeMat;
keeperRecNo = intRecNo(finalKeepers, :);

% make it into a csv file
header = 'varNamesGoesHere';
tableTitle = 'recNo';
tablePrinter(header, keeperRecNo, tableTitle, folderName);




