%findMostRecentLegitRun
function mostRecent = findMostRecentLegitFolder(baseFolder, folderNameStartsWith, mustHaveCSV)

% function call example
%  findMostRecentLegitFolder('C:Dropbox\TBproj\code\outputs\Dec12\GeneXDST_2014', 'c06', 'numQalyPpl_postBurnIn.csv');


%get legit folder names and times
allFileList = dir(baseFolder);

%check if legit
legitMarker = zeros(size(allFileList));
fileDate = zeros(size(allFileList));

for fileNum = 1:size(allFileList,1)
    currentFoldName = allFileList(fileNum).name;
    if strncmp(currentFoldName, folderNameStartsWith,size(folderNameStartsWith,2))
        checker = fullfile(baseFolder,currentFoldName, mustHaveCSV);
        existCode = exist(checker,'file');
        if existCode ~= 0
            legitMarker(fileNum) = 1;  %is a legit folder
            fileDate(fileNum) = datenum(currentFoldName(5:end), 'yyyy-mm-dd_HH-MM-SS') ; %get date
        end
    end
end

if all(legitMarker == 0)
    error('No legit run in %s', baseFolder)
end

%get most recent
[trash,Index] = max(fileDate);
mostRecent = fullfile(baseFolder,allFileList(Index).name);

