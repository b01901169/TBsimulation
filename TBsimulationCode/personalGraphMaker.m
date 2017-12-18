clear
clear all
close all
clc

%makes some tables quickly

masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\apr25_2014\';
folderName = strcat(masterFolder,'paperFigures\quickGraph');
%subFolderList = {'base_fullUptake','GeneXallRNTCP_2_2014_fullUptake','GeneXallRNTCP_3_2014_fullUptake','GeneXallRNTCP_2014_fullUptake','base_fullUptake_post2014','GeneXallRNTCP_2_2014_fullUptake_post2014','GeneXallRNTCP_3_2014_fullUptake_post2014','GeneXallRNTCP_2014_fullUptake_post2014'};

subFolderList = {'base','GeneXallRNTCP_2_2014','GeneXallRNTCP_2014','imperfectBase_mod1','imperfectBase_mod2','imperfectBase_mod3'};
%subFolderList = {'base','GeneXallRNTCP_2_2014','GeneXallRNTCP_2014','GeneX_mod4','GeneX_mod5','GeneX_mod6'};

%subFolderList = {'base','GeneXallRNTCP_2_2014','GeneXallRNTCP_2014','GeneX_mod4','GeneX_mod5'};
subFolderList = {'imperfectBase_mod6','imperfectBase_mod7','imperfectBase_mod11','imperfectBase_mod12'};  %'base','GeneXallRNTCP_2014',

baseFolder = fullfile(masterFolder,subFolderList);

%startStr = {'s', 't','m','e','f'};
startStr = { 'f'};
%startStr = { 's'};

cd(masterFolder);
for i = 1:size(startStr,2)
    %initialize
    diagPoolAll = [];
    numQalysAll = [];
    folderNames = {};
    costsOverTimeAll = [];
    
    for j = 1:size(subFolderList,2)
        cd(masterFolder);        
        folderNameStartsWith = startStr{i};  %this is the first letter of the runs
        mustHaveCSV = 'diagnosedPool.csv';
        folderNameStr = findMostRecentLegitFolder(baseFolder{j}, folderNameStartsWith, mustHaveCSV);
        cd(folderNameStr)
        
        diagPool = dlmread(strcat('diagnosedPool', '.csv'), ',' ,1,0);
        diagPool = diagPool(1560:end,:);
        diagPool = diagPool(:,[2 3 4 6]);
        diagPoolAll = [diagPoolAll, -99901*ones(size(diagPool,1),1), [1996:1/12:2026]' ,diagPool];
        figure(1)
        subplot(1,size(subFolderList,2),j);
        plot([1996:1/12:2026]',diagPool);
        ylabel('diagPool')
        if j == 1
            %legend('activeAll Ppl','coveraledAllPpl_maySeek','coveredAllPpl','canBeDetected','detectedByTest','actuallyDetected');
            legend('coveraledAllPpl_maySeek','coveredAllPpl','canBeDetected','actuallyDetected');
        end
        %ylim([[0 3000]])
        
        
        numQalys = dlmread(strcat('numQalyPpl_postBurnIn', '.csv'), ',' ,1,0);
        numQalys = numQalys(:,1:9);
        numQalysAll = [numQalysAll, -99901*ones(size(numQalys,1),1) ,numQalys];
        figure(2)
        subplot(1,size(subFolderList,2), j);
        plot([1996:1/12:2026]',numQalys(:,1:end));
        ylabel('numQalys')
        if j == 1
            legend('healthyLatNeverTrt','healthyLatPastTrt','healthyLatDOTS','healthyLatMDRTrt','DSpplNoTrt','DSpplInDOTS','DSpplInMDR','MDRpplNotInMDR','MDRpplInMDR');
        end
       
        costsOverTime = dlmread(strcat('costsOverTime', '.csv'), ',' ,1,0);
        costsOverTimeAll = [costsOverTimeAll,costsOverTime];
        figure(3)
        subplot(1,size(subFolderList,2),j);
        plot([1996:1/12:2026]',costsOverTime);
        ylabel('costsOverTime')
        
    end
    
    %print
    cd(masterFolder);
    tableHeader = strcat('diagPoolAll',startStr{i});
    title = strcat('diagPoolAll_',startStr{i});
    tablePrinter(tableHeader, diagPoolAll, title, folderName);
    
    cd(masterFolder);
    tableHeader = strcat('numQalysAll',startStr{i});
    title = strcat('numQalysAll_',startStr{i});
    tablePrinter(tableHeader, numQalysAll, title, folderName);
    
    figure(4)
    titleStr = {'healthyLatNeverTrt','healthyLatPastTrt','healthyLatDOTS','healthyLatMDRTrt','DSpplNoTrt','DSpplInDOTS','DSpplInMDR','MDRpplNotInMDR','MDRpplInMDR'};
    for colInd = 1:9
        colNum = colInd-1;
        %subplot(1,9, colInd);
        figure
        plot([1996:1/12:2026]',numQalysAll(:,[2+colNum: 10 :  (10*size(subFolderList,2))  ]));
        ylabel(titleStr{colInd});
    end
    legStr = {num2str(1)};
    for folds = 2:size(subFolderList,2)
       legStr = {legStr{1:end},num2str(folds)}; 
    end
    legend(legStr)
    
    cd(masterFolder);
    title = strcat('costsOverTimeAll_',startStr{i});
    tableHeader = strcat('costsOverTimeAll', startStr{i});
    tablePrinter(tableHeader, costsOverTimeAll, title, folderName);
    
    i
    figure
    plot([1996:1/12:2026]',costsOverTimeAll)
    ylabel(startStr{i})
    
end


disp('Done!')

