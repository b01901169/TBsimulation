runFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\june11_2014\PPM0p7_0p6_GeneXDST_2014';
masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\june11_2014\';



comp=[];
sensitivityStrList = {'p01','p02','p03'};
for i = 1:3
    cd(masterFolder)
    folderNameStartsWith = sensitivityStrList{i};  %this is the first letter of the runs
    mustHaveCSV = 'monthlyActOutcomes.csv';
    folderNameStr = findMostRecentLegitFolder(runFolder, folderNameStartsWith, mustHaveCSV)
    
    
    cd(folderNameStr)
    monthlyActOutcomes = dlmread(strcat('monthlyActOutcomes', '.csv'), ',' ,1,0);
    comp = [comp,monthlyActOutcomes(:,8)];
    
end
comp