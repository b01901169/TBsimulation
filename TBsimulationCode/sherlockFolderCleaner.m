% sherlockFolderCleaner.m
% remove folders from unfinished runs
clear
clc

masterFold = pwd;
mustHaveCSV = 'monthlyActOutcomes.csv';
TBmac_folders;

for runFoldNum = 1:length(foldNameList)
    
    %go into a folder
    baseFolder = fullfile(masterFold, foldNameList(runFoldNum));
    allFileList = dir(baseFolder{1});
    
    %check if legit
    legitMarker = zeros(size(allFileList));
    fileDate = zeros(size(allFileList));
    
    for runNum = 1:size(allFileList,1)
        currentFoldName = allFileList(runNum).name;
        checker = fullfile(baseFolder,currentFoldName, mustHaveCSV);
        existCode = exist(checker{1},'file');
        if existCode ~= 0
            legitMarker(runNum) = 1;  %is a legit folder
        end
    end
    
    isLegit = logical(legitMarker);
    notLegitIndex = find(~isLegit);
    for notLegitFileNum = 3:size(notLegitIndex,1)  %start at three since the first two are alway . and ..
        thisFile = fullfile(baseFolder,allFileList(notLegitIndex(notLegitFileNum)).name);
        existCode2 = exist(thisFile{1});
        if existCode2 == 7  %delete this file
            fprintf('Deleted: %s \n',thisFile{1})
            rmdir(thisFile{1},'s')            
        end
    end 
end

