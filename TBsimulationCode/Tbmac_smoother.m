%user inputs
minFolders = 3;
maxFolders = 11;
numFolders = maxFolders - minFolders + 1;

masterFold = pwd;
paperFold = fullfile(masterFold, 'paperFigures');

%get the names of folders
TB_MACtargetFolderStrList = {};
index = 1;
for number = minFolders:maxFolders
    TB_MACtargetFolderStrList{index} = strcat('p',sprintf('%02d', number));
    index = index + 1;
end

%go into the correct folder, grab output
for foldNum = 1:length(TB_MACtargetFolderStrList)
    TB_MACtargetFolderStr = TB_MACtargetFolderStrList{foldNum};
    cd(fullfile(paperFold,TB_MACtargetFolderStr));
    calib{foldNum} = dlmread('TB_MAC_calib_results.csv', ',' ,1,0);
    epi{foldNum} = dlmread('TB_MAC_epi_results_all.csv', ',' ,1,3);
    daly1{foldNum} = dlmread('TB_MAC_daly1_results_all.csv', ',' ,1,3);
    daly2{foldNum} = dlmread('TB_MAC_daly2_results_all.csv', ',' ,1,3);
    daly3{foldNum} = dlmread('TB_MAC_daly3_results_all.csv', ',' ,1,4);
    econ{foldNum} = dlmread('TB_MAC_econ_results_all.csv', ',' ,1,3);
    
    warning('off')
    epiHeaders = readtable('TB_MAC_epi_results_all.csv');
    daly1Headers = readtable('TB_MAC_daly1_results_all.csv');
    daly2Headers = readtable('TB_MAC_daly2_results_all.csv');
    daly3Headers = readtable('TB_MAC_daly3_results_all.csv');
    econHeaders = readtable('TB_MAC_econ_results_all.csv');
    warning('on')
    
    foldNum
end

%average across runs (p01-p10)
mean_calib = mean(reshape([calib{:}],size(calib{1},1),size(calib{1},2),numFolders),3);
mean_epi = mean(reshape([epi{:}],size(epi{1},1),size(epi{1},2),numFolders),3);
mean_daly1 = mean(reshape([daly1{:}],size(daly1{1},1),size(daly1{1},2),numFolders),3);
mean_daly2 = mean(reshape([daly2{:}],size(daly2{1},1),size(daly2{1},2),numFolders),3);
mean_daly3 = mean(reshape([daly3{:}],size(daly3{1},1),size(daly3{1},2),numFolders),3);
mean_econ = mean(reshape([econ{:}],size(econ{1},1),size(econ{1},2),numFolders),3);


%smooth over time if needed



%write out
mkdir(fullfile(paperFold,'means'))
cd(fullfile(paperFold,'means'))
writetable(table(mean_calib),'mean_calib.csv');
writetable([epiHeaders(:,1:3),table(mean_epi)],'mean_epi.csv');
writetable([daly1Headers(:,1:3),table(mean_daly1)],'mean_daly1.csv');
writetable([daly2Headers(:,1:3),table(mean_daly2)],'mean_daly2.csv');
writetable([daly3Headers(:,1:4),table(mean_daly3)],'mean_daly3.csv');
writetable([econHeaders(:,1:3),table(mean_econ)],'mean_econ.csv');
cd(masterFold);

%plot
tbmac_plotMaker

