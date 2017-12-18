CEAscenarioFolders = {
    'base','c02';
    'GeneXDST_2014','c04';
    'GeneXallRNTCP_2014','c06';
    'PPM1p0_2014','c08';
    'PPM0p5_2014','c10';
    'PPM1p0_GeneXallRNTCP_2014','c12';
    'PPM0p5_GeneXallRNTCP_2014','c14';
    'PPM0p5_GeneXDST_2014','c16';
    'PPM1p0_GeneXDST_2014','c18';
    };

scriptFolder = pwd;

for scenFoldNum = 1:size(CEAscenarioFolders,1)
    %for scenFoldNum = 1
    baseFolder = fullfile(scriptFolder,CEAscenarioFolders{scenFoldNum,1});
    folderNameStartsWith = CEAscenarioFolders{scenFoldNum,2};
    mustHaveCSV = 'numQalyPpl_postBurnIn.csv';
    mostRecent = findMostRecentLegitFolder(baseFolder, folderNameStartsWith, mustHaveCSV);
    fileToUse{scenFoldNum} = mostRecent;
end

