% sherlockFolderRenamer.m
% rename P02 to p something else
clear
clc

masterFold = pwd;
mustHaveCSV = 'monthlyActOutcomes.csv';
TBmac_folders;
folderNameStartsWith = 'p02';
folderNameChangeTo = 'p11';

for runFoldNum = 1:length(foldNameList)
    
    %go into a folder
    baseFolder = fullfile(masterFold, foldNameList(runFoldNum));
    allFileList = dir(baseFolder{1});
    
    %check if legit
    legitMarker = zeros(size(allFileList));
    fileDate = zeros(size(allFileList));
    
    
    
    for runNum = 1:size(allFileList,1)
        currentFoldName = allFileList(runNum).name;
        
        if strncmp(currentFoldName, folderNameStartsWith,size(folderNameStartsWith,2))
            modFoldName = strcat(folderNameChangeTo, currentFoldName(4:end));
            oldFile = fullfile(baseFolder, currentFoldName)
            newFile = fullfile(baseFolder, modFoldName);            
            movefile(oldFile{:}, newFile{:});            
        end
    end
end

