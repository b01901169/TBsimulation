%rename the file folders

masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\june11_2014\';
subFolderList = {'base','GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};

%get sub folder names
sensStr = {'p','b'};
findFile = {'numQalyPpl_postBurnIn.csv','numQalyPpl_males.csv'};

for subFoldNum = 1:size(subFolderList,2)
    subFolder = fullfile(masterFolder, subFolderList{subFoldNum});
    
    for pbType = 1:2
        cd(masterFolder)
        allLegitRuns = findAllLegitRuns(subFolder, sensStr{pbType} , findFile{pbType});
        folderList = allLegitRuns{2};
        for i = 1:size(folderList,1)
            k = strfind(folderList{i}, '\');
            oldFoldName = folderList{i}((k(end)+1):end);
            newFoldName = strcat( 'p01',oldFoldName(4:end) );
            cd(subFolder);
            if strcmp(oldFoldName, newFoldName) == 0
                movefile(  oldFoldName  ,newFoldName);
            end
        end
    end
    
end