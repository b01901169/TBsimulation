%prevalence by age plot for figure
clear
clear all
printFormat = '-dpdf';
printFormat = '-dpng';

figure1 = 0;
figure2 = 0;
prevIncDoublePlot = 0;  %PAPER FIG 1 prevalence and incidence and MDR plot.  the simulation MDR numbers are read using MDR_whoComparisonPlots.m
makeCasesAverted = 0;  %these appear in the base case output folder as casesAverted2012_2036, deathsAverted, and QALYsGained
makeDemogSummaryTable = 0;  %these appear in the base case output folder as proportionsTables_
burdenFigure = 0;
burdenFigureAgeOnly = 0;
doublePrevGraph = 0;  %PAPER FIG 2
findTotInciCases = 0; %PAPER NUMBER generates the total number of incident cases for the paper
mdrAnnualGraph = 0;
transTxArea = 0; %PAPER FIG 3.  manually put in foldernames.  Smooths over 5 base cases runs.
incidenceByTrans = 0;
organizePlots = 0;  %PAPER FIG 4 this makes the summary subplots graphs when index num is 18.
%need to make foldernames.txt first by using cmd.  cd to folder, then "dir /A:D /S >foldernames.txt", no quotes
%Makes SI figure for half ramp up when index num is 19 and 20
giantPrevGraph = 0;  %%%this does not work yet
ratioOfTransToTx = 0;  %ratio of Trans to Tx graph, gotta manually put in the folder name
totalNumIncidentCases = 0;  %needed for paper..maybe not.  delete this, prolly.
sensi_MDRseed = 0; %SENSITIVITY FIG 1.  must manually put in folder names in section below. what happens with mdr seed different
ageStructure = 0;  %APPENDIX FIGURE: AGE STRUCTURE.  Uses
prevComparisonGraph = 0; %SMDM to compare health outcomes for whatever folders, manual inputs
SMDMpubPrivateGraphs = 0;  %SMDM to make graphs for ppt, manual inputs
prevComparisonTables = 0; %SMDM generates tables for prev and inc
sensitivityAllMDR = 0;  %APPENDIX sensitivity analysis MDR summary graph.  Does not include the vynnicky since it explodes.
make2013PropTransGenMDR = 0;  %this makes the proportion of new cases (incidence) that is transmission generated.  This is a number needed in the paper
TB_MACtargets = 1;  %makes outputs for TB MAC
TB_MACtargets_compiled = 1;  %compiles outputs for TB MAC

%FOLDERNAMES
masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\july8_2014_TBMAC\';  %%MAKE SURE THIS IS CORRECT
%masterFolder = 'C:\Users\ssuen\Desktop\May23p5_2014';

outputFolder = strcat(masterFolder, 'paperFigures');
plotFixerFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\'; %this is where plotfixer lives;

%for FIGURE 1 and FIGURE 2 and appendix
baseStr = 'base\p01_2014-05-28_16-10-36\';
%baseStr = 'GeneXallRNTCP_2014\p01_2014-05-22_13-17-07\';
baseFolder = strcat(masterFolder,baseStr);
noTreatFolder = strcat(masterFolder, 'base\p01_2014-05-28_16-10-36');

%for the TB_MACtargets
TB_MACtargetFolderStr = 'p01';
TB_MACtargetFolds = {
    'base\',...
    'TBMac1_increaseDots_0p25\',...
    'TBMac1_increaseDots_0p5\',...
    'TBMac1_increaseDots_0p75\',...
    'TBMac1_increaseDots\',...
    'TBMac1_increaseDots_adv\',...
    'TBMac2_improveDotsDefault_0p25\',...
    'TBMac2_improveDotsDefault_0p5\',...
    'TBMac2_improveDotsDefault_0p75\',...
    'TBMac2_improveDotsDefault\',...
    'TBMac2_improveDotsDefault_adv\',...
    'TBMac2a_improveDotsDefault_0p25\',...
    'TBMac2a_improveDotsDefault_0p5\',...
    'TBMac2a_improveDotsDefault_0p75\',...
    'TBMac2a_improveDotsDefault\',...
    'TBMac2a_improveDotsDefault_adv\',...
    'TBMac2b_improveDotsDefault_0p25\',...
    'TBMac2b_improveDotsDefault_0p5\',...
    'TBMac2b_improveDotsDefault_0p75\',...
    'TBMac2b_improveDotsDefault\',...
    'TBMac2b_improveDotsDefault_adv\',...
    'TBMac2c_improveDotsDefault_0p25\',...
    'TBMac2c_improveDotsDefault_0p5\',...
    'TBMac2c_improveDotsDefault_0p75\',...
    'TBMac2c_improveDotsDefault\',...
    'TBMac2c_improveDotsDefault_adv\',...
    'TBMac3_XpertReplacesSmear_0p25\',...
    'TBMac3_XpertReplacesSmear_0p5\',...
    'TBMac3_XpertReplacesSmear_0p75\',...
    'TBMac3_XpertReplacesSmear\',...
    'TBMac3_XpertReplacesSmear_adv\',...
    'TBMac3b_XpertReplacesSmear_SSsensUp_0p25\'...
    'TBMac3b_XpertReplacesSmear_SSsensUp_0p5\'...
    'TBMac3b_XpertReplacesSmear_SSsensUp_0p75\'...
    'TBMac3b_XpertReplacesSmear_SSsensUp\'...
    'TBMac3b_XpertReplacesSmear_SSsensUp_adv\'...
    'TBMac4_activeCaseFinding_0p25\',...
    'TBMac4_activeCaseFinding_0p5\',...
    'TBMac4_activeCaseFinding_0p75\',...
    'TBMac4_activeCaseFinding\',...
    'TBMac4_activeCaseFinding_adv\',...
    'TBMac5_ACFandLatentTrt_adv\',...
    'TBMac6_Combination_0p25\',...
    'TBMac6_Combination_0p5\',...
    'TBMac6_Combination_0p75\',...
    'TBMac6_Combination\',...
    'TBMac6_Combination_adv\',...
    };
%     


%TB_MACtargetFolds = {'TBMac3_XpertReplacesSmear','TBMac6_Combination\','TBMac6_combination2\'};
% TB_MACtargetFolds = {...
%     'TBMac6_Combination\',...
% };
%for the MDR calibration plot need to smooth over a few runs
MDRfolders{1,1} = strcat(masterFolder,'base\p01_2014-05-29_18-05-33');
MDRfolders{2,1} = strcat(masterFolder,'base\p01_2014-05-29_18-05-43');
MDRfolders{3,1} = strcat(masterFolder,'base\p01_2014-05-28_16-10-36');
% MDRfolders{4,1} = strcat(masterFolder,'base\s01_2014-04-03_14-40-09');

%for figure3 need many runs of base case to make averages
foldName_paperFig3{1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_1997\a03_2014-01-11_12-28-21';
foldName_paperFig3{2} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2007\a03_2014-01-11_14-20-56';
foldName_paperFig3{3} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2017\a03_2014-01-11_16-13-42';
foldName_paperFig3{4} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2027\a03_2014-01-11_18-06-45';
foldName_paperFig3{5} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2017\a03_2014-01-10_16-32-47';
foldName_paperFig3{6} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2027\a03_2014-01-10_18-29-33';

% foldName_paperFig3{1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_1997\b01_2013-02-13_18-12-58';
% foldName_paperFig3{2} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_2017\b01_2013-02-13_18-21-28';
% foldName_paperFig3{3} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_2027\b01_2013-02-14_00-46-02';
% foldName_paperFig3{4} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\base_60perc_private0p5\b2_2013-02-11_15-21-51';
% foldName_paperFig3{5} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_2007\b01_2013-02-14_13-57-48';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for figure 2
aveTreat_fullCat4_empUptake = strcat(masterFolder, baseStr);
noEvolve = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve\2012-07-02_18-53-03');
aveTreat_noEvolve_2019start = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\June27\inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve2019\2012-07-02_20-47-52';
aveTreat_noEvolve_2027start = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\June27\inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve2027\2012-07-02_22-42-31';
yearPop = [...
    1965    496400000   ;...
    1970    553874000   ;...
    1975    622097000   ;...
    1980    700059000   ;...
    1985    784491000   ;...
    1990    873785000   ;...
    1995    964486000   ;...
    2000    1053898000  ;...
    2005    1140043000  ;...
    2010    1224614000  ;...
    2015    1308221000  ;...
    2020    1386909000  ;...
    2025    1458958000  ;...
    2030    1523482000  ;...
    2035    1579802000  ;...
    2040    1627029000  ;...
    ];
%=======================================================
%=======================================================
%=======================================================
if figure1 == 1
    close all
    %%get data
    cd(baseFolder);
    base_totalPopulationMat = dlmread('HealthOutcomes.csv', ',' ,1,0);  %young
    base{1} = dlmread('actTBoutcomes_age0_20Only.csv', ',' ,1,0);  %young
    base{2} = dlmread('actTBoutcomes_age21_60Only.csv', ',' ,1,0); %mid
    base{3} = dlmread('actTBoutcomes_age61_100Only.csv', ',' ,1,0); % old
    base_smokers{1} = dlmread('actTBoutcomes_smoking_age0_20Only.csv', ',' ,1,0);
    base_smokers{2} = dlmread('actTBoutcomes_smoking_age21_60Only.csv', ',' ,1,0);
    base_smokers{3} = dlmread('actTBoutcomes_smoking_age61_100Only.csv', ',' ,1,0);
    
    cd(noTreatFolder);
    noTreat_PopulationMat = dlmread('HealthOutcomes.csv', ',' ,1,0);  %young
    noTreat{1} = dlmread('actTBoutcomes_age0_20Only.csv', ',' ,1,0);  %young
    noTreat{2} = dlmread('actTBoutcomes_age21_60Only.csv', ',' ,1,0); %mid
    noTreat{3} = dlmread('actTBoutcomes_age61_100Only.csv', ',' ,1,0); % old
    noTreat_smokers{1} = dlmread('actTBoutcomes_smoking_age0_20Only.csv', ',' ,1,0);
    noTreat_smokers{2} = dlmread('actTBoutcomes_smoking_age21_60Only.csv', ',' ,1,0);
    noTreat_smokers{3} = dlmread('actTBoutcomes_smoking_age61_100Only.csv', ',' ,1,0);
    
    %columns are:
    %actSens     actSensCatI     actSensCatII    actSensCatIV    actMdr
    %actMdrCatI  actMdrCatII     actMdrCatIV and then the proportions
    
    %rows are
    %row 7 = 1996, row 10  = 2011, row 15 is 2036.  In five year brackets that
    %end in 1996 (2011, 2036, etc)
    
    startYr = 7;
    endYr = 15;
    
    base_totalPop = sum(base_totalPopulationMat(26:34,1:5),2);
    noTreat_totalPop = sum(noTreat_PopulationMat(26:34,1:5),2);
    
    %get total in each age group
    for i = 1:3
        %loop over young, mid, old to get the total active population in each age
        base_totals{i} = yearPop(startYr:endYr, 2) .* sum(base{i}(startYr:endYr, 1:8),2) ./ base_totalPop;
        noTreat_totals{i} = yearPop(startYr:endYr, 2) .* sum(noTreat{i}(startYr:endYr, 1:8),2) ./ noTreat_totalPop;
    end
    
    %%%%%%%%% treated and untreated %%%%%%%%%%%
    
    %get total treated in each age group
    for i = 1:3
        %loop over young, mid, old to get treated and untreated in every age.
        %i.e.  base_treated{1} is a matrix where rows are 5-years, first col is
        %num treated, and 2nd col is num untreated
        base_treated{i} = [yearPop(startYr:endYr, 2) ,yearPop(startYr:endYr, 2) ] .* [sum(base{i}(startYr:endYr, [2,4,6,7,8]),2) , sum(base{i}(startYr:endYr, [1,5]),2)]./ [base_totalPop,base_totalPop] ;
        noTreat_treated{i} = [yearPop(startYr:endYr, 2) ,yearPop(startYr:endYr, 2) ] .* [sum(noTreat{i}(startYr:endYr, [2,4,6,7,8]),2) , sum(noTreat{i}(startYr:endYr, [1,5]),2) ] ./ [noTreat_totalPop,noTreat_totalPop];
    end
    
    %%%%%%% DS and MDR TB %%%%%%%%%
    
    %get total DS and MDR in each age group
    for i = 1:3
        %loop over young, mid, old to get act TB type in every age.
        %i.e.  base_treated{1} is a matrix where rows are 5-years, first col is
        %num DS, and 2nd col is num MDR
        base_TBtype{i} = [yearPop(startYr:endYr, 2) ,yearPop(startYr:endYr, 2) ] .* [sum(base{i}(startYr:endYr, 1:4),2) , sum(base{i}(startYr:endYr, 5:end),2)] ./ [base_totalPop,base_totalPop];
        noTreat_TBtype{i} = [yearPop(startYr:endYr, 2) ,yearPop(startYr:endYr, 2) ] .* [sum(noTreat{i}(startYr:endYr, 1:4),2) , sum(noTreat{i}(startYr:endYr, 5:end),2) ]./ [noTreat_totalPop,noTreat_totalPop];
    end
    
    %%%%%% smokers and nonsmokers %%%%%%%%%%
    
    %base_totalSmokers should have young matrix in array position {1}, where
    %matrix is rows of 5-years and col1 = num smokers col2 = num nonsmokers
    
    %get total smokers in each age group
    for i = 1:3
        %loop over young, mid, old to get the total active population in each age
        base_totalSmokers{i} = yearPop(startYr:endYr, 2) .* sum(base_smokers{i}(startYr:endYr, 1:8),2)./ [base_totalPop];
        noTreat_totalSmokers{i} = yearPop(startYr:endYr, 2) .* sum(noTreat_smokers{i}(startYr:endYr, 1:8),2)./noTreat_totalPop;
    end
    
    %get total nonsmokers in each age group
    for i = 1:3
        %loop over young, mid, old to get the total active population in each age
        base_totalSmokers{i}(:,2) =  (yearPop(startYr:endYr, 2) .* base_totals{i}./ base_totalPop)   - base_totalSmokers{i} ;
        noTreat_totalSmokers{i}(:,2) = (yearPop(startYr:endYr, 2) .* ( noTreat_totals{i} ./ noTreat_totalPop ) ) - noTreat_totalSmokers{i};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %totals: choices
    %bar([base_totals{1},base_totals{2},base_totals{3}], 'stack')
    %bar([base_treated{1},base_treated{2},base_treated{3}], 'stack')
    %bar([base_treated{1},base_treated{2},base_treated{3}], 'stack')
    %bar([base_totalSmokers{1},base_totalSmokers{2},base_totalSmokers{3}], 'stack')
    
    %NO TREAT (COUNTERFACTUAL)
    bar([1996:5:2036], [noTreat_totals{1},noTreat_totals{2},noTreat_totals{3}], 'stack')
    P=findobj(gca,'type','patch');
    C={[0.6, 0.2, 0.4]  , [0.9, 0.8, 0.4], [0.5, 0.7, 0.9]  }; % make a colors list
    for numI = 1:length(P)
        set(P(numI),'facecolor',C{numI});
    end
    ylim([0 12000000])
    %legend('threslow','threshigh','vratio','Location','NorthWest')
    pdfTitle = 'noTreat_actTB';
    print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
    
    bar([1996:5:2036], [base_totals{1},base_totals{2},base_totals{3}], 'stack')
    P=findobj(gca,'type','patch');
    C={[0.6, 0.2, 0.4], [0.9, 0.8, 0.4], [0.5, 0.7, 0.9]   }; % make a colors list
    for j = 1:length(P)
        set(P(j),'facecolor',C{j});
    end
    ylim([0 12000000])
    pdfTitle = 'base_actTB';
    print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
    
    close all
    base_allAges = base_totals{1} + base_totals{2} + base_totals{3};
    noTreat_allAges = noTreat_totals{1} + noTreat_totals{2} + noTreat_totals{3};
    plot([1996:5:2036],base_allAges, 'Color', [0.2, 0.3, 1])
    hold on;
    plot([1996:5:2036],noTreat_allAges, 'Color', [0.8, 0.2, 0.4])
    legend('with DOTS','without DOTS','Location','NorthWest')
    pdfTitle = 'base_actTB_line';
    print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
    
    close all
    base_allAges = base_totals{1} + base_totals{2} + base_totals{3};
    noTreat_allAges = noTreat_totals{1} + noTreat_totals{2} + noTreat_totals{3};
    plot([1996:5:2036],noTreat_allAges, 'Color', [0.8, 0.2, 0.4])
    hold on;
    area([1996:5:2036], base_totals{1}+base_totals{2}+base_totals{3}, 'Facecolor', [0.5, 0.7, 0.9])
    hold on;
    area([1996:5:2036], base_totals{1}+base_totals{2}, 'Facecolor', [0.2, 0.5, 1])
    hold on;
    area([1996:5:2036], base_totals{1}, 'Facecolor', [0.2, 0.3, 1])
    legend('without DOTS, all ages','with DOTS, active TB age 60+', 'with DOTS, active TB age 21-60', 'with DOTS, active TB age 0-20','Location','NorthWest')
    pdfTitle = 'base_actTB_line_Age';
    print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
    
    close all
    %cases averted by age
    averted{1} = noTreat_totals{1} - base_totals{1};
    averted{2} = noTreat_totals{2} - base_totals{2};
    averted{3} = noTreat_totals{3} - base_totals{3};
    area([1996:5:2036], averted{1}+averted{2}+averted{3}, 'Facecolor', [0.5, 0.7, 0.9])
    hold on;
    area([1996:5:2036], averted{1}+averted{2}, 'Facecolor', [0.2, 0.5, 1])
    hold on;
    area([1996:5:2036], averted{1}, 'Facecolor', [0.2, 0.3, 1])
    ylim([0 8000000])
    legend('averted active TB age 60+', 'averted active TB age 21-60', 'averted active TB age 0-20','Location','NorthWest')
    pdfTitle = 'base_actTB_avertedByAge';
    print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
    
    
    %     interleaved = reshape([noTreat_allAges;base_allAges],9,[]);
    %     h = bar([1996:5:2036],interleaved);
    %     set(get(h(1),'BaseLine'));
    %     P=findobj(gca,'type','patch');
    %     C={[0.5, 0.7, 0.9], [0.9, 0.8, 0.4] }; % make a colors list
    %     legend('no DOTS', 'with DOTS', 'location', 'NorthWest');
    %     pdfTitle = 'interleavedBars';
    %     print('-dpdf','-r100',[outputFolder '/' pdfTitle]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if figure2 == 1
    %%%FIGURE 2
    % clearvars -except aveTreat_fullCat4_empUptake aveTreat_noEvolve_2019start aveTreat_noEvolve_2027start noEvolve outputFolder
    close all
    
    startYear = 26;
    
    cd(aveTreat_fullCat4_empUptake);
    aveTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    baseCase_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
    cd(aveTreat_noEvolve_2019start)
    aveTreat_noEvolve_2019start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    noEvolve_2019start_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
    cd(aveTreat_noEvolve_2027start)
    aveTreat_noEvolve_2027start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    noEvolve_2027start_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
    cd(noEvolve)
    noEvolve_2011start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    noEvolve_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
    
    
    
    %%Graph for noEvolve starting in 2011, 2019, 2027
    for typeTB = 4  %loop over latDS, actDS, latMDR, actMDR
        
        if typeTB == 1
            healthOutcomesCol = 7; %do latent DS
            TBtypeTitleStr = 'Latent DS';
            ylimMat = [0 0.8];
        elseif typeTB == 2
            healthOutcomesCol = 9; %active DS
            TBtypeTitleStr = 'Active DS';
            ylimMat = [0 0.012];
        elseif typeTB == 3
            healthOutcomesCol = 8; %latent MDR
            TBtypeTitleStr = 'Latent MDR';
            ylimMat = [0 0.01];
        elseif typeTB == 4
            healthOutcomesCol = 10; %act MDR
            TBtypeTitleStr = 'Active MDR';
            ylimMat = [0 0.0002];
        end
        
        base = aveTreat_fullCat4_empUptake_outcomes(startYear:34, healthOutcomesCol);
        start2011 =  noEvolve_2011start_outcomes(startYear:34, healthOutcomesCol);
        start2019 = aveTreat_noEvolve_2019start_outcomes(startYear:34, healthOutcomesCol);
        start2027 = aveTreat_noEvolve_2027start_outcomes(startYear:34, healthOutcomesCol);
        
        years = [1996:5:2036]';
        plot(years,  smooth(start2027,3,'moving') , 'Color', [0.6 0 0.6], 'LineWidth', 2 );
        hold on;
        plot(years, smooth(start2019,3,'moving'), 'Color', [0.9 0.4 0.2], 'LineWidth', 2);
        hold on;
        plot(years, smooth(start2011,3,'moving') , 'Color',  [0.9 0.7 0.1], 'LineWidth', 2);
        hold on;
        plot(years, smooth(base,3,'moving'), 'Color',  [0.1 0.5 0.9], 'LineWidth', 2);
        legend('No Treatment-generated MDR Starting 2027', 'No Treatment-generated MDR Starting 2019','No Treatment-generated MDR Starting 2011','No Intervention', 'location', 'NorthWest')
        ylim([0 0.0003])
        %    cd(plotFixerFold);
        %    plotfixer;
        printStr = strcat('delayInNoEvolve_', TBtypeTitleStr);
        print('-dpdf','-r100',[outputFolder '/' printStr]);
        
    end
    
    
    %%%%%%graph 2a
    close all
    months = [1996:1/12:(2036+11/12)]';
    
    extraOverlayGraph = 0;
    yearly_MDRevolvedTrans = zeros(floor(((2036-1996+1)*12)/60),2);  %because numNewMDR_evolvedOrTrans_postBurnIn is in months from 1996 to 2046, so need to stop at 492, the end of 2036
    %sum to get five year bins
    for i = 0:(floor(((2036-1996+1)*12)/60))
        yearly_MDRevolvedTrans(i+1,:) = sum( baseCase_MDRevolvedTrans( (60*i+1):60*i+60, 1:2) );
    end
    
    totNumPpl = sum(aveTreat_fullCat4_empUptake_outcomes(startYear:33,1:6)')';
    totNewMDRcases = sum(yearly_MDRevolvedTrans')';
    percEvolved = yearly_MDRevolvedTrans(:,1)./totNewMDRcases(:,1);
    percTransmitted = yearly_MDRevolvedTrans(:,2)./totNewMDRcases(:,1);
    
    toPlot = [yearly_MDRevolvedTrans(:,1)./totNumPpl, yearly_MDRevolvedTrans(:,2)./totNumPpl];
    
    plot(years(1:size(years,1),:), smooth(percTransmitted, 3, 'moving'), 'LineWidth', 2,'Color', [0.9 0.4 0.2] );
    hold on;
    plot(years(1:size(years,1),:), smooth(percEvolved, 3, 'moving'), 'LineWidth', 2,'Color',  [0.6 0 0.6]);
    legend('Transmission Generated', 'Treatment Generated', 'location', 'NorthWest')
    % xlim([2008 2035]);
    ylim([0 1]);
    fileStr = strcat('MDRevolvedTransmitted');
    print('-dpdf','-r100',[outputFolder '/' fileStr]);
    
    close all
    India2036pop = yearPop(end,2);
    bar(India2036pop.*[ sum(baseCase_MDRevolvedTrans(481, 1:2))/baseCase_MDRevolvedTrans(481,3);...
        sum(noEvolve_MDRevolvedTrans(481, 1:2))/noEvolve_MDRevolvedTrans(481,3);...
        sum(noEvolve_2019start_MDRevolvedTrans(481, 1:2))/noEvolve_2019start_MDRevolvedTrans(481,3);...
        sum(noEvolve_2027start_MDRevolvedTrans(481, 1:2))/noEvolve_2027start_MDRevolvedTrans(481,3) ]);
    fileStr = strcat('MDRevolvedTransmitted_overTime');
    print('-dpdf','-r100',[outputFolder '/' fileStr]);
    
    close all
    plotyears = interp1(yearPop(:,1), yearPop(:,2), [2011:1/12:2036] );
    base =  repmat(plotyears, 2,1)'.*baseCase_MDRevolvedTrans(181:481, 1:2)./baseCase_MDRevolvedTrans(481,3);
    noEvolve_2011 = repmat(plotyears, 2, 1)'.*noEvolve_MDRevolvedTrans(181:481, 1:2)./noEvolve_MDRevolvedTrans(481,3);
    noEvolve_2019 = repmat(plotyears, 2, 1)'.*noEvolve_2019start_MDRevolvedTrans(181:481, 1:2)./noEvolve_2019start_MDRevolvedTrans(481,3);
    noEvolve_2027 = repmat(plotyears, 2, 1)'.*noEvolve_2027start_MDRevolvedTrans(181:481, 1:2)./noEvolve_2027start_MDRevolvedTrans(481,3);
    averted2011 = sum(base, 1) - sum(noEvolve_2011, 1) ;
    averted2019 = sum(base, 1) - sum(noEvolve_2019, 1) ;
    averted2027 = sum(base, 1) - sum(noEvolve_2027, 1) ;
    
    bar([2011, 2019, 2027], [averted2011; averted2019 ; averted2027])
    P=findobj(gca,'type','patch');
    C={[0.6, 0.2, 0.4], [0.5, 0.7, 0.9]   }; % make a colors list
    for j = 1:length(P)
        set(P(j),'facecolor',C{j});
    end
    fileStr = 'MDRaverted';
    legend('averted treatment-generated cases', 'averted transmitted cases')
    print('-dpdf','-r100',[outputFolder '/' fileStr]);
    
    
end
%=======================================================
%=======================================================
%=======================================================
if prevIncDoublePlot == 1
    
    currentDirectory = pwd;
    cd(baseFolder);
    plotResolution = '-r100';
    
    %INCIDENCE
    
    %who incidence data
    dataTime = [1990    ; 1991  ; 1992  ; 1993  ; 1994  ;  1995 ; 1996  ; 1997  ; 1998  ; 1999  ; 2000  ; 2001  ; 2002  ; 2003  ; 2004  ; 2005  ; 2006  ; 2007  ; 2008  ; 2009  ; 2010  ];
    actTBinci = [216    ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 216   ; 215   ; 214   ; 212   ; 209   ; 205   ; 201   ; 196   ; 190   ; 185 ];
    distToUpperLimit = [39  ; 37    ; 35    ; 33    ; 32    ; 30    ; 29    ; 27    ; 26    ; 25    ; 24    ; 23    ; 23    ; 23    ; 22    ; 22    ; 22    ; 21    ; 21    ; 21    ; 20    ];
    distToLowerLimit = [35  ;  33   ; 32    ; 30    ; 29    ; 27    ; 26    ; 25    ; 24    ; 23    ; 22    ; 22    ; 22    ; 22    ; 22    ;21     ; 21    ; 20    ; 20    ; 19    ; 18    ];
    %simulation data
    yrlyPlottingMat = dlmread('TBincidence.csv', ',' ,1,0);
    toPlot = yrlyPlottingMat(120:end,1:2) ./ [yrlyPlottingMat(120:end,3) ,yrlyPlottingMat(120:end,3)  ] ;
    %make legend
    smallx = [1996:1:2011];
    toPlotLeg = toPlot * 100000;
    p = plot(smallx,toPlotLeg(11:26));
    set(p,'Color','red','LineWidth',1);
    hold on;
    errorbar(dataTime(7:end), actTBinci(7:end),distToLowerLimit(7:end),distToUpperLimit(7:end));
    hleg1 = legend('Simulation Estimates','WHO Estimates','Orientation','horizontal'); %'Location','South'
    set(hleg1, 'Position', [.35,.15,.1,.2])
    cd(plotFixerFold);
    plotFixerPNAS;
    cd(currentDirectory);
    print('-dpdf',plotResolution,[outputFolder '/' 'calibrationTriplePanel_legend.pdf']);
    close all
    subplot(2,3,[2 5]);
    %make simulation incidences into rates out of 100000 and plot
    smallx = [1996:1:2011];
    toPlot = toPlot * 100000;
    %p = plot(x,toPlot);
    p = plot(smallx,toPlot(11:26));
    set(p,'Color','red','LineWidth',1);
    %title({'Simulation Comparisons With WHO estimates','Cases Out of 100000 Individuals'})
    xlabel('Year');
    ylabel('Active TB Annual Incidence (out of 100000 People)');
    hold on;
    errorbar(dataTime(7:end), actTBinci(7:end),distToLowerLimit(7:end),distToUpperLimit(7:end));
    ylim([0 400]);
    xlim([1995 2012]);
    cd(plotFixerFold);
    plotFixerPNAS;
    cd(currentDirectory);
    
    %PREVALENCE
    subplot(2,3,[1 4]);
    
    
    
    %who prevalence data
    dataYear = [1990    ; 1991  ; 1992  ; 1993  ; 1994  ;  1995 ; 1996  ; 1997  ; 1998  ; 1999  ; 2000  ; 2001  ; 2002  ; 2003  ; 2004  ; 2005  ; 2006  ; 2007  ; 2008  ; 2009  ; 2010  ];
    actTBPrev = [459 ; 460 ; 460    ; 461   ; 462   ; 462   ; 463   ; 464   ; 464   ; 465   ; 466   ; 466   ; 436   ; 409   ; 383   ; 358   ; 335   ; 314   ; 294   ; 275   ; 256];
    distToUpperLimit = [56  ; 56    ; 56    ; 56    ; 56    ; 57    ; 56    ; 56    ; 57    ; 56    ; 56    ; 57    ; 62    ; 66    ; 71    ; 78    ; 84    ; 91    ; 99    ; 107   ; 117];
    distToLowerLimit = [52  ; 53    ; 52    ; 53    ; 53    ; 53    ; 53    ; 53    ; 53    ; 53    ; 54    ; 53    ; 57    ; 62    ; 66    ; 70    ; 74    ; 80    ; 85    ; 90    ; 95];
    %simulation data
    cd(baseFolder);
    simOut = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
    cd(currentDirectory); %get back out to the original directory
    scaled = simOut./repmat(sum(simOut(:,1:11),2),1,size(simOut,2));
    toPlot = sum(scaled(:,4:11),2);
    %plot
    smallx = [1996:1/12:2010];
    toPlot = (toPlot(1:169,1))*100000;
    p = plot(smallx,toPlot);
    set(p,'Color','red','LineWidth',1);
    %title({'Prevalence With WHO estimates','Cases Out of 100000 Individuals'})
    xlabel('Year');
    ylabel({'','Active TB Prevalence (out of 100000 People)'});
    hold on;
    errorbar(dataYear(7:end), actTBPrev(7:end),distToLowerLimit(7:end),distToUpperLimit(7:end));
    ylim([0 600]);
    xlim([1995 2012]);
    cd(plotFixerFold);
    plotFixerPNAS;
    cd(currentDirectory);
    
    %MDR graph
    simMDRoutput = MDR_whoComparisonPlots(MDRfolders);
    
    subplot(2,3,6);
    %these numbers are read from an excel sheet in Dropbox/TBproject/MDRmatching_validation.xslx
    who_percMDRamongNewCases2008 = [2.3 1.8 2.8];
    
    who_percMDRamongNewCases2008 = repmat(who_percMDRamongNewCases2008,3,1); %gotta hack it since the horizontal bars are weird and can't edit the errorbar.m file
    %OLD  sim_calc_percMDRinNewCases2008 = [1.46    1.30    1.61];  %this is where the MDR_whoComparisonPlots outputs are manually put in
    % sim_calc_percMDRinNewCases2008 = [1.3 1.1 1.5];  %base case
    sim_calc_percMDRinNewCases2008 = simMDRoutput(1,1:3)*100;
    %sim_calc_percMDRinNewCases2008 = [1.06441084798753       0.958557410124268         1.1702642858508];  %scenario for minor revision plos one
    sim_calc_percMDRinNewCases2008 = repmat(sim_calc_percMDRinNewCases2008, 3,1);
    
    %plot
    errorbar([2006.9 2007.9 2009.1], who_percMDRamongNewCases2008(:,1),who_percMDRamongNewCases2008(:,1)-who_percMDRamongNewCases2008(:,2),who_percMDRamongNewCases2008(:,3)-who_percMDRamongNewCases2008(:,1),'bo');
    hold on;
    errorbar([2006.9 2008.1 2009.1] , sim_calc_percMDRinNewCases2008(:,1), sim_calc_percMDRinNewCases2008(:,1)- sim_calc_percMDRinNewCases2008(:,2),sim_calc_percMDRinNewCases2008(:,3)-sim_calc_percMDRinNewCases2008(:,1),'ro');
    xlim([2007 2009]);
    ylim([0 5]);
    set(gca, 'XTick', [2008]);
    ylabel('% MDR among new TB cases');
    xlabel('Year');
    cd(plotFixerFold);
    plotFixerPNAS;
    cd(currentDirectory);
    
    subplot(2,3,3);
    who_percMDRamongNewCases2008=[ 99000   79000   120000];
    who_percMDRamongNewCases2008 = repmat(who_percMDRamongNewCases2008,3,1); %gotta hack it since the horizontal bars are weird and can't edit the errorbar.m file
    %OLD  sim_calc_percMDRinNewCases2008=[  136561  84613   188509];
    %sim_calc_percMDRinNewCases2008=[ 96415.6901375837          77234.5021520272           115596.87812314];  %base case
    sim_calc_percMDRinNewCases2008 = simMDRoutput(4,1:3);
    %sim_calc_percMDRinNewCases2008=[71593.9190403257          67954.9628086836          75232.8752719678];  %scenario for minor revision plos one
    
    sim_calc_percMDRinNewCases2008 = repmat(sim_calc_percMDRinNewCases2008, 3,1);
    errorbar([2006.9 2007.9 2009.1], who_percMDRamongNewCases2008(:,1),who_percMDRamongNewCases2008(:,1)-who_percMDRamongNewCases2008(:,2),who_percMDRamongNewCases2008(:,3)-who_percMDRamongNewCases2008(:,1),'bo');
    hold on;
    errorbar([2006.9 2008.1 2009.1], sim_calc_percMDRinNewCases2008(:,1), sim_calc_percMDRinNewCases2008(:,1)- sim_calc_percMDRinNewCases2008(:,2),sim_calc_percMDRinNewCases2008(:,3)-sim_calc_percMDRinNewCases2008(:,1),'ro');
    xlim([2007 2009]);
    ylim([0 250000]);
    set(gca, 'XTick', [2008]);
    ylabel({'No. Incident MDR-TB Cases'});
    % set(gcf, 'PaperPositionMode', 'manual');
    % set(gcf, 'PaperUnits', 'inches');
    % set(gcf, 'PaperPosition', [0.1 0.1 3.42 3.42]);
    cd(plotFixerFold);
    plotFixerPNAS;
    cd(currentDirectory);
    print(printFormat,plotResolution,[outputFolder '/' 'calibrationTriplePanel']);
    close all
    
    %     set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]); %default for PaperPosition
    %     set(gcf, 'PaperUnits', 'inches');  %set back to default
    
    cd(currentDirectory);
end
%=======================================================
%=======================================================
%=======================================================
if makeCasesAverted == 1
    cd(masterFolder);
    casesAvertedMaker(masterFolder, 1,noTreatStr, baseStr, 'compare to no treat')
end
%=======================================================
%=======================================================
%=======================================================
if makeDemogSummaryTable ==1
    
    
    for type = 1:2
        currentFold = pwd;
        cd(baseFolder);
        if type == 1
            %male female
            type1_age{1} = dlmread('HealthOutcomes_postBurnIn_male_age0_20Only.csv', ',' ,1,0);%young
            type1_age{2} = dlmread('HealthOutcomes_postBurnIn_male_age21_60Only.csv', ',' ,1,0);%mid
            type1_age{3} = dlmread('HealthOutcomes_postBurnIn_male_age61_100Only.csv', ',' ,1,0); %old
            type2_age{1} = dlmread('HealthOutcomes_postBurnIn_female_age0_20Only.csv', ',' ,1,0);  %young
            type2_age{2} = dlmread('HealthOutcomes_postBurnIn_female_age21_60Only.csv', ',' ,1,0); %mid
            type2_age{3} = dlmread('HealthOutcomes_postBurnIn_female_age61_100Only.csv', ',' ,1,0);  %old
            typeStr = 'sex';
            header = 'young_actTB, young_actMDR, mid_actTB, mid_actMDR, old_actTB, old_actMDR, col are men2011 men2036 women2011 women2036';
        elseif type == 2
            %smoker nonsmoker
            type1_age{1} = dlmread('HealthOutcomes_postBurnIn_smoking_age0_20Only.csv', ',' ,1,0);%young
            type1_age{2} = dlmread('HealthOutcomes_postBurnIn_smoking_age21_60Only.csv', ',' ,1,0);%mid
            type1_age{3} = dlmread('HealthOutcomes_postBurnIn_smoking_age61_100Only.csv', ',' ,1,0); %old
            type2_age{1} = dlmread('HealthOutcomes_postBurnIn_nonSmoking_age0_20Only.csv', ',' ,1,0);  %young
            type2_age{2} = dlmread('HealthOutcomes_postBurnIn_nonSmoking_age21_60Only.csv', ',' ,1,0); %mid
            type2_age{3} = dlmread('HealthOutcomes_postBurnIn_nonSmoking_age61_100Only.csv', ',' ,1,0);  %old
            typeStr = 'smokStat';
            header = 'young_actTB, young_actMDR, mid_actTB, mid_actMDR, old_actTB, old_actMDR, col are smoker2011 smoker2036 nonsmoker2011 nonsmoker2036';
        end
        cd(currentFold);
        
        for i = 1:3
            table(1,1+2*(i-1)) = type1_age{i}(10,9); %prop actDS TB 2011
            table(1,2+2*(i-1)) = type1_age{i}(10, 10); %prop actMDR TB 2011
            table(2,1+2*(i-1)) = type1_age{i}(15,9); %prop actDS TB 2036
            table(2,2+2*(i-1)) = type1_age{i}(15, 10); %prop actMDR TB 2036
        end
        
        for i = 1:3
            table(3,1+2*(i-1)) = type2_age{i}(10,9); %prop actDS TB 2011
            table(3,2+2*(i-1)) = type2_age{i}(10, 10); %prop actMDR TB 2011
            table(4,1+2*(i-1)) = type2_age{i}(15,9); %prop actDS TB 2036
            table(4,2+2*(i-1)) = type2_age{i}(15, 10); %prop actMDR TB 2036
        end
        name = strcat('proportionsTables_', typeStr);
        tablePrinter(header, table, name, baseFolder);
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make a burden figure
if burdenFigure == 1
    currentFold = pwd;
    cd(baseFolder);
    TBtypeNames = {'latSens','latMDR','actSens','actMDR'};
    all_out =  dlmread('HealthOutcomes.csv', ',' ,1,0); %everyone
    for type = 1:2
        if type == 1
            cd(baseFolder);
            %male female
            type1_age{1} = dlmread('HealthOutcomes_male_age0_20Only.csv', ',' ,1,0);%young
            type1_age{2} = dlmread('HealthOutcomes_male_age21_60Only.csv', ',' ,1,0);%mid
            type1_age{3} = dlmread('HealthOutcomes_male_age61_100Only.csv', ',' ,1,0); %old
            type2_age{1} = dlmread('HealthOutcomes_female_age0_20Only.csv', ',' ,1,0);  %young
            type2_age{2} = dlmread('HealthOutcomes_female_age21_60Only.csv', ',' ,1,0); %mid
            type2_age{3} = dlmread('HealthOutcomes_female_age61_100Only.csv', ',' ,1,0);  %old
            typeDesc = 'sex';
        elseif type == 2
            cd(baseFolder);
            %smoker nonsmoker
            type1_age{1} = dlmread('HealthOutcomes_smoking_age0_20Only.csv', ',' ,1,0);%young
            type1_age{2} = dlmread('HealthOutcomes_smoking_age21_60Only.csv', ',' ,1,0);%mid
            type1_age{3} = dlmread('HealthOutcomes_smoking_age61_100Only.csv', ',' ,1,0); %old
            type2_age{1} = dlmread('HealthOutcomes_nonSmoking_age0_20Only.csv', ',' ,1,0);  %young
            type2_age{2} = dlmread('HealthOutcomes_nonSmoking_age21_60Only.csv', ',' ,1,0); %mid
            type2_age{3} = dlmread('HealthOutcomes_nonSmoking_age61_100Only.csv', ',' ,1,0);  %old
            typeDesc = 'smokeStat';
        end
        cd(currentFold);
        
        for TBtype = 2:5
            plot(:,1) = type1_age{3}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            plot(:,2) = type1_age{2}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            plot(:,3) = type1_age{1}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            plot(:,4) = type2_age{1}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            plot(:,5) = type2_age{2}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            plot(:,6) = type2_age{3}(34:-1:29,TBtype)./all_out(34:-1:29, TBtype);
            typeStr = TBtypeNames{TBtype-1};
            
            barh([2036:-5:2011],plot, 'stack')
            xlim([0,1]);
            if type == 1
                leg = legend('male 60+','male 21-60', 'male 20-', 'female 20-', 'female 21-60', 'female 60+','Location','NorthEastOutside');
            elseif type == 2
                leg = legend('smoker 60+','smoker 21-60', 'smoker 20-', 'nonsmoker 20-', 'nonsmoker 21-60', 'nonsmoker 60+','Location','NorthEastOutside');
            end
            set(leg);
            fileNum = TBtype-1;
            fileStr = strcat('burdenPlot_',typeStr, '_',typeDesc);
            print('-dpdf','-r100',[outputFolder '/' fileStr]);
        end
    end
    
end
if burdenFigureAgeOnly == 1
    currentFold = pwd;
    cd(baseFolder);
    TBtypeNames = {'latSens','latMDR','actSens','actMDR'};
    all_out =  dlmread('HealthOutcomes.csv', ',' ,1,0); %everyone
    type1_age{1} = dlmread('HealthOutcomes_age0_20Only.csv', ',' ,1,0);%young
    type1_age{2} = dlmread('HealthOutcomes_age21_60Only.csv', ',' ,1,0);%mid
    type1_age{3} = dlmread('HealthOutcomes_age61_100Only.csv', ',' ,1,0); %old
    cd(currentFold);
    for TBtype = 2:5
        plot(:,1) = type1_age{1}(29:34,TBtype)./all_out(29:34, TBtype);
        plot(:,2) = type1_age{2}(29:34,TBtype)./all_out(29:34, TBtype);
        plot(:,3) = type1_age{3}(29:34,TBtype)./all_out(29:34, TBtype);
        typeStr = TBtypeNames{TBtype-1};
        barh([2011:5:2036],plot, 'stack')
        xlim([0,1]);
        set(gca, 'YDir', 'reverse');
        %         P=findobj(gca,'type','patch');
        %         C={[0.6, 0.2, 0.4], [0.9, 0.8, 0.4], [0.5, 0.7, 0.9]   }; % make a colors list
        %         for j = 1:length(P)
        %             set(P(j),'facecolor',C{j});
        %         end
        leg = legend('Aged 20-','Aged 21-60', 'Aged 60+','Location','NorthEastOutside');
        set(leg);
        fileNum = TBtype-1;
        fileStr = strcat('AgeBurdenPlot_',typeStr);
        print('-dpdf','-r100',[outputFolder '/' fileStr]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make Jeremy's prev fig2
if doublePrevGraph == 1
    currentFold = pwd;
    cd(baseFolder);
    aveTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    cd(noTreatFolder);
    noTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
    cd(currentFold);
    
    plotResolution = '-r100';
    
    for startYrLoop = 2:2
        years = [2011:5:2036]';
        startYear = 29;
        printStr = 'start2011_2013';
        startMon = 205;
        startCalYear = 2013;
        
        if startYrLoop == 2
            years = [1996:5:2036]';
            startYear = 26;
            printStr = 'start1996';
            startMon = 1;
            startCalYear = 1996;
        end
        
        %get data for plotting
        DSpropAveTreat = aveTreat_fullCat4_empUptake_outcomes(startYear:34, 9);  %act DS
        DSpropNoTreat = noTreat_fullCat4_empUptake_outcomes(startYear:34, 9);  %act DS
        MDRpropAveTreat = aveTreat_fullCat4_empUptake_outcomes(startYear:34, 10); %actMDR
        MDRpropNoTreat = noTreat_fullCat4_empUptake_outcomes(startYear:34, 10); %actMDR
        
        
        %         subplot(2,1,1);
        %         [AX,H1,H2] = plotyy(years,DSpropNoTreat,years,MDRpropNoTreat,'plot');
        %         set(get(AX(1),'Ylabel'),'String','Non-MDR Prevalence (% of total pop)');
        %         set(get(AX(2),'Ylabel'),'String','MDR Prevalence (% of total pop)');
        %         set(H1,'LineStyle','--');
        %         set(H2,'LineStyle','--');
        %         set(AX(1),'ylim',[0 0.008],'ytick',[0:0.002:0.008]);
        %         set(AX(2),'ylim',[0 0.0004],'ytick',[0:0.0001:0.0004]);
        %
        %         hold on;
        %         [BX,J1,J2] = plotyy(years,DSpropAveTreat,years,MDRpropAveTreat,'plot');
        %         set(J1,'LineStyle','-');
        %         set(J2,'LineStyle','-');
        %         set(BX(1),'ylim',[0 0.008],'ytick',[0:0.002:0.008]);
        %         set(BX(2),'ylim',[0 0.0004],'ytick',[0:0.0001:0.0004]);
        %         title('Disease Prevalence With and Without DOTS');
        %make the legend
        plot(years, DSpropNoTreat, 'k--');
        hold on;
        plot(years, DSpropAveTreat, 'k');
        hold on;
        plot([2013 2013], [0 0.005], 'k:');
        ylim([0 0.009]);
        hleg1=legend('Without Treatment','With Treatment','Current Time: 2013','Orientation', 'horizontal'); % 'Location','NorthWest',
        set(hleg1, 'Position', [.45,.75,.1,.2])
        cd(plotFixerFold);
        plotFixerPNAS;
        cd(currentFold);
        print('-dpdf',plotResolution,[outputFolder '/' 'prevAverDS_MDR_legend']);
        close all
        
        %Color for light blue [0 0.75 0.95],  dark blue [0.1 0.1 0.75]
        subplot(2,2,1);
        plot(years, DSpropNoTreat, '--', 'Color', [0 0.75 0.95]);
        hold on;
        plot(years, DSpropAveTreat, 'Color', [0 0.75 0.95]);
        hold on;
        plot([2013 2013], [0 0.008], 'k:');
        ylabel('Non-MDR Prevalence'); % (% of Total Pop)
        xlabel('Year');
        title('Non-MDR');
        xlim([1996 2038]);
        cd(plotFixerFold);
        plotFixerPNAS;
        cd(currentFold);
        
        subplot(2,2,2);
        plot(years, MDRpropNoTreat, '--','Color', [0.85 0.3 0.4]);
        hold on;
        plot(years, MDRpropAveTreat, 'Color', [0.85 0.3 0.4]);
        hold on;
        plot([2013 2013], [0 0.0006], 'k:');
        ylabel('MDR Prevalence');  % (% of Total Pop)
        xlabel('Year');
        title('MDR');
        ylim([0 0.0006]);
        xlim([1996 2038]);
        cd(plotFixerFold);
        plotFixerPNAS;
        cd(currentFold);
        
        currentFold = pwd;
        cd(baseFolder);
        aveTreat_deaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
        cd(noTreatFolder);
        noTreat_deaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
        cd(currentFold);
        
        %start from 1 if starting jan 1996
        %181 if jan 2011, 492 for dec 2036 ; 205 if jan 2013, 516 for dec 2038
        endMon = 516;
        endCalYear = 2038;
        
        %Load historic population in India
        IndianPop = [...
            1950    371857000   ;...
            1955    406374000   ;...
            1960    447844000   ;...
            1965    496400000   ;...
            1970    553874000   ;...
            1975    622097000   ;...
            1980    700059000   ;...
            1985    784491000   ;...
            1990    873785000   ;...
            1995    964486000   ;...
            2000    1053898000  ;...
            2005    1140043000  ;...
            2010    1224614000  ;...
            2015    1308221000  ;...
            2020    1386909000  ;...
            2025    1458958000  ;...
            2030    1523482000  ;...
            2035    1579802000  ;...
            2040    1627029000  ;...
            2045    1664519000  ;...
            2050    1692008000  ;...
            ];
        shortIndianPop = interp1(IndianPop(:,1), IndianPop(:,2), [startCalYear:1/12:endCalYear+(11/12)])';
        
        %make the cumulative deaths from disease over time
        aveTreat_cumuDeaths = zeros(endMon - startMon+1, 2);
        noTreat_cumuDeaths = zeros(endMon - startMon+1, 2);
        
        aveTreat_runningTot = aveTreat_deaths(startMon, 1:2)/aveTreat_deaths(startMon, 3)*shortIndianPop(1);
        aveTreat_cumuDeaths(1,1:2)  = aveTreat_runningTot;
        noTreat_runningTot = noTreat_deaths(startMon, 1:2)/noTreat_deaths(startMon, 3)*shortIndianPop(1);
        noTreat_cumuDeaths(1,1:2) = noTreat_runningTot ;
        for monLoop = startMon+1 :endMon
            aveTreat_cumuDeaths(monLoop - startMon +1, :) = aveTreat_runningTot + (aveTreat_deaths(monLoop, 1:2)./ repmat(aveTreat_deaths(monLoop,3),1,2) .* repmat(shortIndianPop(monLoop-startMon+1),1,2) );
            aveTreat_runningTot = aveTreat_cumuDeaths(monLoop - startMon+1, :);
            noTreat_cumuDeaths(monLoop - startMon+1, :) = noTreat_runningTot + (noTreat_deaths(monLoop, 1:2)./ repmat(noTreat_deaths(monLoop,3),1,2) .* repmat(shortIndianPop(monLoop-startMon+1),1,2) );
            noTreat_runningTot = noTreat_cumuDeaths(monLoop - startMon+1, :);
        end
        avertedDeaths = noTreat_cumuDeaths - aveTreat_cumuDeaths;
        
        %plot
        toPlot = -avertedDeaths;
        %   toPlot(:,2) = -toPlot(:,2);
        years = [startCalYear:1/12:endCalYear+(11/12)];
        
        %         subplot(2,1,2);
        %         plot(years', [toPlot(:,1)/1000000,toPlot(:,2)/1000000]);
        %         hold on;
        %         plot(years, [toPlot(:,2)/1000000], 'Color', [0.3, 0.7, 0.5]);
        %         legend('DOTS averted non-MDR deaths', 'Additional MDR deaths with DOTS', 'Location', 'NorthWest');
        %         ylim([0 70]);
        %         title(sprintf('Cumulative Deaths Since %i', startCalYear));
        %         ylabel('In Millions');
        %         xlabel('Year');
        % %AVERTED VERSION
        %         toPlot = -avertedDeaths;
        %         years = [startCalYear:1/12:endCalYear+(11/12)];
        %         subplot(2,2,3);
        %         plot(years, 20*ones(size(years)), 'k--');
        %         hold on;
        %         plot(years, toPlot(:,1)/1000000, 'k');
        %         hold on;
        %         plot([2013 2013], [-100 20], 'k:');
        %         legend('Without DOTS','With DOTS','Current Time: 2013', 'Location','SouthWest');
        %
        %         subplot(2,2,3);
        %         plot(years, toPlot(:,1)/1000000, 'Color', [0 0.75 0.95]);
        %         ylabel('Number of Deaths, In Millions');
        %         hold on;
        %         plot([2013 2013], [-100 20], 'k:');
        %         xlabel('Year');
        %         title('Additional non-MDR deaths with DOTS');
        %
        %
        %         subplot(2,2,4);
        %         plot(years, toPlot(:,2)/1000000, 'Color', [0.85 0.3 0.4]);
        %         hold on;
        %         plot([2013 2013], [0 3], 'k:');
        %         ylim([0 3]);
        %         ylabel('Number of Deaths, In Millions');
        %         xlabel('Year');
        %         title('Additional MDR deaths with DOTS');
        subplot(2,2,3);
        plot(years, noTreat_cumuDeaths(:,1)/1000000, '--','Color', [0 0.75 0.95]);
        hold on;
        plot(years, aveTreat_cumuDeaths(:,1)/1000000,'Color', [0 0.75 0.95]);
        ylabel({'Deaths Since 1996 (Millions)'});
        hold on;
        plot([2013 2013], [-100 150], 'k:');
        xlabel('Year');
        %title('Non-MDR deaths');
        ylim([0 150]);
        xlim([1996 2038]);
        cd(plotFixerFold);
        plotFixerPNAS;
        cd(currentFold);
        
        subplot(2,2,4);
        plot(years, noTreat_cumuDeaths(:,2)/1000000, '--','Color', [0.85 0.3 0.4]);
        hold on;
        plot(years, aveTreat_cumuDeaths(:,2)/1000000,'Color', [0.85 0.3 0.4]);
        hold on;
        plot([2013 2013], [0 6], 'k:');
        ylim([0 6]);
        ylabel({'Deaths Since 1996 (Millions)'});
        xlabel('Year');
        %title('MDR deaths');
        xlim([1996 2038]);
        
        %save the graph
        %   cd(plotFixerFold);
        %   plotfixer;
        plotName = strcat('prevAverDS_MDR_', printStr);
        cd(plotFixerFold);
        plotFixerPNAS;
        cd(currentFold);
        print(printFormat,plotResolution,[outputFolder '/' plotName]);
        close all
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%transmission vs evolved graph
if mdrAnnualGraph == 1
    for dotsPlusOrNot = 1:4
        %%find right folders to read from
        if dotsPlusOrNot == 1
            noEvolveFolds = cell(27,1);
            noEvolveFolds{13} = '2012-07-02_18-53-03';
            noEvolveFolds{14} = '2012-07-03_00-37-20';
            noEvolveFolds{15} = '2012-07-03_02-32-09';
            noEvolveFolds{16} = '2012-07-03_04-26-57';
            noEvolveFolds{17} = '2012-07-03_06-21-48';
            noEvolveFolds{18} = '2012-07-03_08-16-01';
            noEvolveFolds{19} = '2012-07-02_20-47-52';
            noEvolveFolds{20} = '2012-07-03_12-20-06';
            noEvolveFolds{21} = '2012-07-03_14-16-06';
            noEvolveFolds{22} = '2012-07-03_16-11-59';
            noEvolveFolds{23} = '2012-07-03_18-07-57';
            noEvolveFolds{24} = '2012-07-03_20-03-38';
            noEvolveFolds{25} = '2012-07-03_21-55-28';
            noEvolveFolds{26} = '2012-07-03_23-46-35';
            noEvolveFolds{27} = '2012-07-02_22-42-31';
            dotsStr = 'wDotsPlus';
        elseif dotsPlusOrNot == 2
            noEvolveNoIVFolds{13} = '2012-07-04_18-00-38';
            noEvolveNoIVFolds{14} = '2012-07-04_14-33-36';
            noEvolveNoIVFolds{15} = '2012-07-02_18-22-28';
            noEvolveNoIVFolds{16} = '2012-07-03_16-20-43';
            noEvolveNoIVFolds{17} = '2012-07-03_16-17-40';
            noEvolveNoIVFolds{18} = '2012-07-02_18-18-28';
            noEvolveNoIVFolds{19} = '2012-07-03_10-24-13';
            noEvolveNoIVFolds{20} = '2012-07-02_18-17-45';
            noEvolveNoIVFolds{21} = '2012-07-02_18-16-34';
            noEvolveNoIVFolds{22} = '2012-07-02_18-15-32';
            noEvolveNoIVFolds{23} = '2012-07-03_20-15-34';
            noEvolveNoIVFolds{24} = '2012-07-03_20-14-30';
            noEvolveNoIVFolds{25} = '2012-07-03_20-13-46';
            noEvolveNoIVFolds{26} = '2012-07-04_19-51-46';
            noEvolveNoIVFolds{27} = '2012-07-03_20-11-30';
            dotsStr = 'noDotsPlus';
        elseif dotsPlusOrNot == 3
            noEvolveNoIVALLFolds{13} = '2012-07-04_21-43-27';
            noEvolveNoIVALLFolds{14} = '2012-07-04_23-34-38';
            noEvolveNoIVALLFolds{15} = '2012-07-05_01-26-03';
            noEvolveNoIVALLFolds{16} = '2012-07-05_03-17-06';
            noEvolveNoIVALLFolds{17} = '2012-07-05_05-08-43';
            noEvolveNoIVALLFolds{18} = '2012-07-05_07-00-26';
            noEvolveNoIVALLFolds{19} = '2012-07-04_18-19-04';
            noEvolveNoIVALLFolds{20} = '2012-07-05_08-51-49';
            noEvolveNoIVALLFolds{21} = '2012-07-04_18-15-46';
            noEvolveNoIVALLFolds{22} = '2012-07-04_18-14-47';
            noEvolveNoIVALLFolds{23} = '2012-07-04_18-12-28';
            noEvolveNoIVALLFolds{24} = '2012-07-04_18-12-19';
            noEvolveNoIVALLFolds{25} = '2012-07-04_18-10-54';
            noEvolveNoIVALLFolds{26} = '2012-07-04_18-10-11';
            noEvolveNoIVALLFolds{27} = '2012-07-04_18-08-24';
            dotsStr = 'noDotsPlusUniverse';
        elseif dotsPlusOrNot == 4
            noEvolvehalfIVFolds{13} = '2012-07-06_01-00-15';
            noEvolvehalfIVFolds{14} = '2012-07-06_00-07-52';
            noEvolvehalfIVFolds{15} = '2012-07-05_23-44-51';
            noEvolvehalfIVFolds{16} = '2012-07-05_23-58-17';
            noEvolvehalfIVFolds{17} = '2012-07-05_23-47-03';
            noEvolvehalfIVFolds{18} = '2012-07-05_19-49-19';
            noEvolvehalfIVFolds{19} = '2012-07-05_23-28-34';
            noEvolvehalfIVFolds{20} = '2012-07-05_19-47-18';
            noEvolvehalfIVFolds{21} = '2012-07-05_19-45-05';
            noEvolvehalfIVFolds{22} = '2012-07-05_19-43-42';
            noEvolvehalfIVFolds{23} = '2012-07-05_19-42-11';
            noEvolvehalfIVFolds{24} = '2012-07-05_19-37-10';
            noEvolvehalfIVFolds{25} = '2012-07-05_19-36-49';
            noEvolvehalfIVFolds{26} = '2012-07-05_19-34-18';
            noEvolvehalfIVFolds{27} = '2012-07-05_19-41-23';
            dotsStr = 'halfDotsPlus';
        end
        %%%READ IN DATA
        currentFold = pwd;
        if dotsPlusOrNot <= 3
            cd(baseFolder);
        elseif dotsPlusOrNot == 4
            cd(strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta\', '2012-07-05_19-22-25'));
        end
        
        aveTreatMDR = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        %read in all the noEvolve data
        noEvolve = cell(27,1);
        if dotsPlusOrNot == 1
            cd(strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve\',noEvolveFolds{13}));      %for 2013 noEvolve
        elseif dotsPlusOrNot == 2
            cd(strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_noCatIV_empUpta_noEvolve2013\',noEvolveNoIVFolds{13}));      %for 2013 no DOTS+, noEvolve
        elseif dotsPlusOrNot == 3
            cd(strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_ALLnoCatIV_empUpta_noEvolve2013\',noEvolveNoIVALLFolds{13}));      %for 2013 no DOTS+ ever, noEvolve
        elseif dotsPlusOrNot == 4
            cd(strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta_noEvolve2013\',noEvolvehalfIVFolds{13}));      %for 2013 no DOTS+ ever, noEvolve
        end
        
        noEvolve{13} = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        for i = 14:27  %for 2014 to 2027 noEvolve
            if dotsPlusOrNot == 1
                cd(strcat(masterFolder,sprintf('inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve20%i', i),'\', noEvolveFolds{i} ));
            elseif dotsPlusOrNot == 2
                cd(strcat(masterFolder,sprintf('inf0p0023_lat2p16_aveTreat_noCatIV_empUpta_noEvolve20%i', i),'\', noEvolveNoIVFolds{i} ));      %for 2013 no DOTS+, noEvolve
            elseif dotsPlusOrNot == 3
                cd(strcat(masterFolder,sprintf('inf0p0023_lat2p16_aveTreat_ALLnoCatIV_empUpta_noEvolve20%i', i),'\', noEvolveNoIVALLFolds{i} ));      %for 2013 no DOTS+, noEvolve
            elseif dotsPlusOrNot == 4
                cd(strcat(masterFolder,sprintf('inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta_noEvolve20%i', i),'\', noEvolvehalfIVFolds{i} ));      %for 2013 no DOTS+, noEvolve
            end
            noEvolve{i} = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        end
        cd(currentFold);
        
        %%SET UP PLOT
        plotResolution = '-r100';
        startCalYear = 2013;
        endCalYear = 2038;
        
        %Load historic population in India
        IndianPop = [...
            1950    371857000   ;...
            1955    406374000   ;...
            1960    447844000   ;...
            1965    496400000   ;...
            1970    553874000   ;...
            1975    622097000   ;...
            1980    700059000   ;...
            1985    784491000   ;...
            1990    873785000   ;...
            1995    964486000   ;...
            2000    1053898000  ;...
            2005    1140043000  ;...
            2010    1224614000  ;...
            2015    1308221000  ;...
            2020    1386909000  ;...
            2025    1458958000  ;...
            2030    1523482000  ;...
            2035    1579802000  ;...
            2040    1627029000  ;...
            2045    1664519000  ;...
            2050    1692008000  ;...
            ];
        shortIndianPop = interp1(IndianPop(:,1), IndianPop(:,2), [startCalYear:1/12:endCalYear+(11/12)])';
        
        %make the cumulative cases of MDR over time
        %initalize
        noEvolve_cumuCases = cell(27,1);
        aveTreat_cumuCases = cell(27,1);
        
        for i = 13:27
            %initalize
            startMon = 205 + (12*(i-13));  %startMon = 205;  %2013
            %endMon = 516;  %endMon = 516;  %2038
            endMon = startMon+120;  %endMon = 516;  %2038
            
            aveTreat_cumuCases{i} = zeros(endMon - startMon+1, 2);
            noEvolve_cumuCases{i} = zeros(endMon - startMon+1, 2);
            
            %make the cumulative total for baseCase
            aveTreat_runningTot = aveTreatMDR(startMon, 1:2)/aveTreatMDR(startMon, 3)*shortIndianPop(1);
            aveTreat_cumuCases{i}(1,1:2)  = aveTreat_runningTot;
            for monLoop = startMon+1 :endMon
                aveTreat_cumuCases{i}(monLoop - startMon +1, :) = aveTreat_runningTot + (aveTreatMDR(monLoop, 1:2)/aveTreatMDR(monLoop, 3)*shortIndianPop(monLoop-startMon+1));
                aveTreat_runningTot = aveTreat_cumuCases{i}(monLoop - startMon+1, :);
            end
            
            %now all the noEvolve scenarios make cumulative total
            noEvolve_runningTot = noEvolve{i}(startMon, 1:2)/noEvolve{i}(startMon, 3)*shortIndianPop(1);
            noEvolve_cumuCases{i}(1,1:2) = noEvolve_runningTot;
            for monLoop = startMon+1 :endMon
                noEvolve_cumuCases{i}(monLoop - startMon+1, :) = noEvolve_runningTot + (noEvolve{i}(monLoop, 1:2)/noEvolve{i}(monLoop, 3)*shortIndianPop(monLoop-startMon+1));
                noEvolve_runningTot = noEvolve_cumuCases{i}(monLoop - startMon+1, :);
            end
        end
        
        %%COMPARE AVE TREAT AND NOEVOLVE SCENARIOS
        toPlot = zeros(27-13,2);
        for i = 13:27
            noEvolve_avertedIndia{i} = aveTreat_cumuCases{i}(end,:) - noEvolve_cumuCases{i}(end,:);
            %toPlot(i-12,:) = [noEvolve_avertedIndia{i}(1), noEvolve_avertedIndia{i}(2)]./repmat((2038-2000-i),1,2);  %if end time point is 2038 for all scenarios
            toPlot(i-12,:) = [noEvolve_avertedIndia{i}(1), noEvolve_avertedIndia{i}(2)]./repmat(10,1,2);    %if end time point is 10 years out for all scenarios
            %toPlot(i-12,:) = [noEvolve_avertedIndia{i}(1), noEvolve_avertedIndia{i}(2)];  %not annualized
        end
        
        %plot
        years = [startCalYear:1:2027]';
        
        %%ONE VERSION
        for versionType = 1:6
            if versionType == 1
                %bar, stacked
                clear numI lengthP P C
                bar(years, toPlot, 'stack')
                P=findobj(gca,'type','patch');
                C={[0.6, 0.2, 0.4], [0.5, 0.7, 0.9] }; % make a colors list
                for numI = 1:length(P)
                    set(P(numI),'facecolor',C{numI});
                end
                ylim([0 500000]);
                versionTypeStr = 'stackedBar';
            elseif versionType == 2
                %bar, nonstacked
                clear numI lengthP P C
                bar(years, [toPlot(:,2),toPlot(:,1)])
                P=findobj(gca,'type','patch');
                C={[0.6, 0.2, 0.4], [0.5, 0.7, 0.9] }; % make a colors list
                for numI = 1:length(P)
                    set(P(numI),'facecolor',C{numI});
                end
                ylim([0 300000]);
                versionTypeStr = 'nonStackedBar';
            elseif versionType == 3
                %smoothed area
                clear numI lengthP P C
                area(years, smooth(toPlot(:,1)+toPlot(:,2),3,'moving'), 'Facecolor', [0.6, 0.2, 0.4])
                hold on;
                area(years,smooth(toPlot(:,1),3,'moving'), 'Facecolor', [0.5, 0.7, 0.9])
                ylim([0 500000]);
                versionTypeStr = 'smoothedArea';
            elseif versionType == 4
                %%NONsmoothed area
                clear numI lengthP P C
                area(years, toPlot(:,1)+toPlot(:,2), 'Facecolor', [0.6, 0.2, 0.4])
                hold on;
                area(years,toPlot(:,1), 'Facecolor', [0.5, 0.7, 0.9])
                ylim([0 500000]);
                versionTypeStr = 'nonsmoothedArea';
            elseif versionType == 5
                %smoothed line
                clear numI lengthP P C
                plot(years, [smooth(toPlot(:,1),10,'moving')/1000, smooth(toPlot(:,2),10,'moving')/1000 ] );
                P=findobj(gca,'type','patch');
                C={[0.6, 0.2, 0.4], [0.5, 0.7, 0.9] }; % make a colors list
                for numI = 1:length(P)
                    set(P(numI),'facecolor',C{numI});
                end
                %ylim([0 300]);
                versionTypeStr = 'smoothedLine';
            elseif versionType == 6
                %dots
                clear numI lengthP P C
                %plot(years, [toPlot(:,1)/1000, toPlot(:,2)/1000 ] );
                %hold on;
                plot(years, [toPlot(:,1)/1000, toPlot(:,2)/1000 ], 'o' );
                ylim([0 300]);
                versionTypeStr = 'unsmoothedLine';
            end
            
            xlabel('Year DOTS improved');
            ylabel('MDR Cases Averted Annually');
            %legend('Averted Treatment-Generated Cases','Averted Transmission-Generated Cases');
            legend('Averted Transmission-Generated Cases', 'Averted Treatment-Generated Cases');
            
            %save the graph
            %   cd(plotFixerFold);
            %   plotfixer;
            typeStr = strcat('evolTrans_MDRaverted_', versionTypeStr,'_', dotsStr);
            print('-dpdf',plotResolution,[outputFolder '/' typeStr]);
            close all
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if transTxArea == 1
    plotResolution = '-r100';
    %cd('aveTreat_aveCatIV_empUpta_slowDST_2007\2012-07-20_18-10-53');
    %cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\July19\aveTreat_aveCatIV_empUpta_slowDST_2007\2012-07-20_18-10-53');
    cd(baseFolder);
    xvals = [1996:1/12:2038+11/12] ;
    
    foldName = foldName_paperFig3;
    
    for fold = 1:size(foldName,2)
        cd(foldName{fold})
        graphDataRaw = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        if fold == 1
            cola = graphDataRaw(:,1);
            colb = graphDataRaw(:,2);
            colc = graphDataRaw(:,3);
            cold = graphDataRaw(:,4);
            cole = graphDataRaw(:,1) + graphDataRaw(:,4);
        elseif fold ~= 1
            cola = [cola, graphDataRaw(:,1)];
            colb = [colb, graphDataRaw(:,2)];
            colc = [colc, graphDataRaw(:,3)];
            cold = [cold, graphDataRaw(:,4)];
            cole = [cola, graphDataRaw(:,1)] + [cold, graphDataRaw(:,4)];
        end
    end
    
    %   graphData = [mean(cola,2),mean(colb,2),mean(colc,2),mean(cold,2)];  %use this if you want 3 colors (DOTS evolved, trans, and privClinic evolved graph).  Uncomment col0 line below.
    graphData = [mean(cole,2),mean(colb,2),mean(colc,2)];  %use this if you only want the actual paper figure (evolved or trans graph)
    %   graphData = [mean(cola,2),mean(colb,2),mean(colc,2)];  %ignoring private clinics.
    
    smoothFactor = 240;  %usually 240
    col0 = zeros(516,1);
    %   col0 = smooth(graphData(1:516,4)./graphData(1:516,3),smoothFactor,'moving');    % keep if you want 3 colors.  Default commented out.
    col2 = smooth(graphData(1:516,2)./graphData(1:516,3),smoothFactor,'moving');
    col3 = smooth(graphData(1:516,1)./graphData(1:516,3),smoothFactor,'moving');
    toPlot = [col0, col0+col3, col0+col3+col2];
    toPlot = toPlot * 12 * 100000;  %if we want annual incidence out of 100000
    
    H0 = area(xvals(1:end-12), toPlot(1:end-12,3 ));
    set(H0,'FaceColor',[0.2 0.7 0.9]);
    hold on;
    
    H1 = area(xvals(1:end-12), toPlot(1:end-12,2 ));
    set(H1,'FaceColor',[1 0.83 0.3]);
    hold on;
    if size(graphData,2) == 4
        H2 = area(xvals(1:end-12), toPlot(1:end-12,1));
        set(H2,'FaceColor',[1 0 0.3]);
    end
    xlabel('Year');
    ylabel('Incidence Rate (annual cases per 100,000)');
    plot([2013, 2013], [0 18], 'k:')
    %plot([2013, 2013], [0 0.000014], 'k:')
    %plot([2013, 2013], [0 0.000002], 'k:')    %smaller y-scale for runs where MDR incidence is low
    legend('Transmission-generated MDR Incidence Rate','Treatment-generated MDR Incidence Rate','Current Time: 2013','Location','NorthWest');
    xlim([1995 xvals(end-12)]);
    
    currentFold = pwd;
    cd(plotFixerFold);
    % plotfixer;
    typeStr = strcat('tran_txGenInci');
    subOutFold = strcat(outputFolder);
    cd(plotFixerFold);
    plotFixerPNAS;
    print(printFormat,plotResolution,[subOutFold '/' typeStr]);
    close all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if incidenceByTrans == 1
    plotResolution = '-r100';
    cd(masterFolder);
    
    %use this if you want all the folders in a directory to be done
    string = textread('folderNames.txt', '%s');
    j = 1;
    for i = 1:size(string,1)
        if strcmp(string{i},'Directory')
            if strcmp(string{i+1},'of')
                if strcmp(string{i+2}(end-2), '-')
                    graphFolder{j} = string(i+2);
                    holdThis = graphFolder{j};
                    startStr{j} = holdThis{1}(end-23:end-20);
                    slashPositions = strfind(holdThis{1}, '\');
                    caseStr{j} = holdThis{1}(slashPositions(1,end-1)+1: slashPositions(1,end)-5 );
                    discontinuity{j} = str2num(startStr{j});
                    if strcmp(caseStr{j}, 'aveTreat_aveCatIV_empUpta_slowDST_')
                        discontinuity{j} = 2038;
                    end
                    j = j+1;
                end
            end
        end
    end
    
    
    %MAKE THE INCIDENCE PLOT
    for j = 1:size(graphFolder,2)
        folderNa = graphFolder{j};
        folderNameStr = folderNa{1}(1:end);
        cd(folderNameStr);
        
        for smoothFlag = 1:2
            if smoothFlag == 1
                smoothFactor = 240; smoothStr = 'smoothed';
            elseif smoothFlag == 2
                smoothFactor = 1; smoothStr = 'unsmoothed';
            end
            
            for MDRtype = 1:2
                if MDRtype == 1
                    graphData = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
                    col1 = (graphData(1:516,1)+graphData(1:516,2))./graphData(1:516,3);
                    col2 = graphData(1:516,2)./graphData(1:516,3);
                    col3 = graphData(1:516,1)./graphData(1:516,3);
                    endTime = 516;
                    breakTime = ((discontinuity{j}-1996)*12)+12;
                    mdrStr = 'MDRincidence';
                    xvals = [1996:1/12:2038+11/12] ;
                    ylimitMat = [0 0.00001];
                elseif MDRtype == 2
                    nonMDRinci = dlmread('TBincidence.csv', ',' ,1,0);
                    graphData = [nonMDRinci(130:172,1), zeros(43,1), nonMDRinci(130:172,3)];
                    col3 = graphData(:,1)./graphData(:,3);
                    col2 = zeros(size(col3));
                    col1 = zeros(size(col3));
                    endTime = size(col1,1);
                    breakTime = (discontinuity{j}-1996);
                    mdrStr = 'nonMDRincidence';
                    xvals = [1996:1:2038];
                    ylimitMat = [0 0.0025];
                end
                
                if smoothFactor ~= 1 && discontinuity{j} == 2038
                    col1 = smooth((graphData(1:endTime,1)+graphData(1:endTime,2))./graphData(1:endTime,3), smoothFactor, 'moving');
                    col2 = smooth(graphData(1:endTime,2)./graphData(1:endTime,3), smoothFactor, 'moving');
                    col3 = smooth(graphData(1:endTime,1)./graphData(1:endTime,3), smoothFactor, 'moving');
                elseif smoothFactor ~= 1
                    timeP = {1, breakTime, breakTime, breakTime, breakTime+1, endTime};
                    
                    col1first = smooth((graphData(timeP{1}:timeP{2},1)+graphData(timeP{1}:timeP{2},2))./graphData(timeP{1}:timeP{2},3), smoothFactor, 'moving');
                    col1sec = []; % (graphData(timeP{3}:timeP{4},1)+graphData(timeP{3}:timeP{4},2))./graphData(timeP{3}:timeP{4},3);
                    col1third = smooth((graphData(timeP{5}:timeP{6},1)+graphData(timeP{5}:timeP{6},2))./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                    col1 = [col1first;col1sec;col1third];
                    
                    col2first = smooth(graphData(1:timeP{2},2)./graphData(1:timeP{2},3), smoothFactor, 'moving');
                    col2sec  =[]; % graphData(timeP{3}:timeP{4},2)./graphData(timeP{3}:timeP{4},3);
                    col2third = smooth(graphData(timeP{5}:timeP{6},2)./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                    col2 = [col2first;col2sec; col2third];
                    
                    col3first = smooth(graphData(1:timeP{2},1)./graphData(1:timeP{2},3), smoothFactor, 'moving');
                    col3sec = []; %graphData(timeP{3}:timeP{4},1)./graphData(timeP{3}:timeP{4},3);
                    col3third = smooth(graphData(timeP{5}:timeP{6},1)./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                    col3 = [col3first;col3sec;col3third];
                end
                
                toPlot = [col1, col2, col3];
                plot(xvals(1:end-24), toPlot(1:end-24,:));
                %         h = findobj(gca,'Type','line');
                %         set(h(1),'Color','k');
                %         set(h(2),'Color',[0.8, 0.2, 0.4]);
                %         set(h(3),'Color',[0.5, 0.6, 1]);
                %
                ylabel({'Incidence Rate', '(%age of population with active Disease)'});
                xlabel('Year');
                if MDRtype == 1
                    legend('Total','Transmission-generated', 'Treatment-generated', 'Location', 'NorthWest');
                elseif MDRtype == 2
                    title(num2str(toPlot(1,end)));
                end
                ylim(ylimitMat);
                
                %save the graph
                currentFold = pwd;
                cd(plotFixerFold);
                plotfixer;
                typeStr = strcat(caseStr{j},startStr{j},'_',smoothStr);
                subOutFold = strcat(outputFolder,'\figure4\',mdrStr,'\',smoothStr,'\', startStr{j});
                print('-dpdf',plotResolution,[subOutFold '/' typeStr]);
                cd(currentFold);
                close all
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if organizePlots == 1
    display('in organizePlots')
    
    plotResolution = '-r300';
    cd(masterFolder);
    
    %initialize graphFolder to have all the plots
    string = textread('folderNames.txt', '%s');
    j = 1;
    for i = 1:size(string,1)
        if strcmp(string{i},'Directory')
            if strcmp(string{i+1},'of')
                if strcmp(string{i+2}(end-2), '-')
                    graphFolder{j} = string(i+2);
                    holdThis = graphFolder{j};
                    startStr{j} = holdThis{1}(end-27:end-24);
                    slashPositions = strfind(holdThis{1}, '\');
                    caseStr{j} = holdThis{1}(slashPositions(1,end-1)+1: slashPositions(1,end)-5 );
                    discontinuity{j} = str2num(startStr{j});
                    j = j+1;
                end
            end
        end
    end
    
    smoothFactor = 24;
    
    %initialize plotNameArray
    %DOTS x DOTS+
    plotNameArray{1} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_perfCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','bestTreat_perfCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_perfCatIV_empUpta_slowDST_noEvolve','','','',6};
    totPlotName{1} = 'DOTS_DOTSplus';
    sizeMat{1,1} = 3; sizeMat{1,2}= 3;
    positionMat{1} = [1,2,4,5,7,8];
    %DOTSx DST
    plotNameArray{2} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve',9};
    totPlotName{2} = 'DOTS_DST';
    sizeMat{2,1} = 3; sizeMat{2,2}= 3;
    positionMat{2} = [1,2,3,4,5,6,7,8,9];
    %DOTSx speed
    plotNameArray{3} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_fastUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_fastUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_fastUpta_slowDST_noEvolve','','','',6};
    totPlotName{3} = 'DOTS_Uptake';
    sizeMat{3,1} = 3; sizeMat{3,2}= 3;
    positionMat{3} = [1,2,4,5,7,8];
    %DOTS+ x DST
    plotNameArray{4} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_perfCatIV_empUpta_slowDST_','aveTreat_perfCatIV_empUpta_fastDST_','aveTreat_perfCatIV_empUpta_accDST','','','',6};
    totPlotName{4} = 'DOTSplus_DST';
    sizeMat{4,1} = 3; sizeMat{4,2}= 3;
    positionMat{4} = [1,2,3,4,5,6];
    %DOTS+ x speed
    plotNameArray{5} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_fastUpta_slowDST_','aveTreat_perfCatIV_empUpta_slowDST_','aveTreat_perfCatIV_fastUpta_slowDST_','','','','','',4};
    totPlotName{5} = 'DOTSplus_speed';
    sizeMat{5,1} = 3; sizeMat{5,2}= 3;
    positionMat{5} = [1,2,4,5];
    %DST x speed
    plotNameArray{6} = {'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_fastUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_fastUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_fastUpta_accDST','','','',6};
    totPlotName{6} = 'DST_Uptake';
    sizeMat{6,1} = 3; sizeMat{6,2}= 3;
    positionMat{6} = [1,2,4,5,7,8];
    %main effects
    plotNameArray{7} = {'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_perfCatIV_empUpta_slowDST_','aveTreat_aveCatIV_fastUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST',7};
    totPlotName{7} = 'mainEffects';
    sizeMat{7,1} = 1; sizeMat{7,2}= 7;
    positionMat{7} = [1,2,3,4,5,6,7];
    %secondaryEffects
    plotIndex = 8;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_fastUpta_fastDST_','aveTreat_aveCatIV_fastUpta_accDST','aveTreat_perfCatIV_empUpta_fastDST_','aveTreat_perfCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_fastUpta_slowDST_noEvolve','aveTreat_perfCatIV_empUpta_slowDST_noEvolve','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_fastUpta_slowDST_','bestTreat_perfCatIV_empUpta_slowDST_',12};
    totPlotName{plotIndex} = 'secondaryEffects';
    sizeMat{plotIndex,1} = 3; sizeMat{plotIndex,2}= 4;
    positionMat{plotIndex} = [1:1:12];
    %tertiaryEffects
    plotIndex = 9;
    plotNameArray{plotIndex}={'aveTreat_perfCatIV_fastUpta_fastDST_','aveTreat_perfCatIV_fastUpta_accDST','aveTreat_aveCatIV_fastUpta_fastDST_noEvolve','aveTreat_aveCatIV_fastUpta_accDST_noEvolve','aveTreat_perfCatIV_empUpta_fastDST_noEvolve','aveTreat_perfCatIV_empUpta_accDST_noEvolve','aveTreat_perfCatIV_fastUpta_slowDST_noEvolve','bestTreat_aveCatIV_fastUpta_fastDST_','bestTreat_aveCatIV_fastUpta_accDST','bestTreat_perfCatIV_empUpta_fastDST_','bestTreat_perfCatIV_empUpta_accDST','bestTreat_perfCatIV_fastUpta_slowDST_',12};
    totPlotName{plotIndex} = 'tertiaryEffects';
    sizeMat{plotIndex,1} = 3; sizeMat{plotIndex,2}= 5;
    positionMat{plotIndex} = [1,2,6,7,8,9,10,11,12,13,14,15];
    %summaryGraph
    plotIndex = 10;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST',20};
    totPlotName{plotIndex} = 'summaryGraphs';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 5;
    positionMat{plotIndex} = [1:1:20];
    sumPlot = 1;
    sumYrChange{plotIndex} = [5:5:20];
    
    %summaryGraph2
    plotIndex = 11;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve',36};
    totPlotName{plotIndex} = 'summaryGraphs2';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 9;
    positionMat{plotIndex} = [1:1:36];
    sumPlot = 1;
    sumYrChange{plotIndex} = [9:9:36];
    
    %1997, 2017, 2027
    plotIndex = 12;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve',54};
    totPlotName{plotIndex} = 'summaryGraphs3';
    sizeMat{plotIndex,1} = 6; sizeMat{plotIndex,2}= 9;
    positionMat{plotIndex} = [1:1:54];
    sumPlot = 1;
    
    sumYrChange{plotIndex} = [9:9:54];
    
    %1997, 2007, 2017, 2027
    plotIndex = 13;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve',72};
    totPlotName{plotIndex} = 'summaryGraphs4';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 9;
    positionMat{plotIndex} = [1:1:72];
    sumPlot = 1;
    sumYrChange{plotIndex} = [9:9:72];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %%%paper figure attempt
    %Main Effects best and noEvolve only 1997, 2007, 2017, 2027
    plotIndex = 14;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_accDST',40};
    totPlotName{plotIndex} = 'summaryMain';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 5;
    positionMat{plotIndex} = [1:1:40];
    sumPlot = 1;
    sumYrChange{plotIndex} = [5:5:40];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %Main Effects best and noEvolve only 1997, 2007, 2017, 2027
    plotIndex = 15;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_accDST',24};
    totPlotName{plotIndex} = 'summaryIntaxn_best';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %Main Effects best and noEvolve only 1997, 2007, 2017, 2027
    plotIndex = 16;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_noEvolve','aveTreat_aveCatIV_empUpta_accDST_noEvolve',24};    totPlotName{plotIndex} = 'summaryIntaxn_noEvolve';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %Main Effects best and noEvolve only 1997, 2007, 2017, 2027
    plotIndex = 17;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_noEvolve','aveTreat_aveCatIV_empUpta_fastDST_',32};
    totPlotName{plotIndex} = 'summarySimpleAllPrev';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 4;
    positionMat{plotIndex} = [1:1:32];
    sumPlot = 1;
    sumYrChange{plotIndex} = [4:4:32];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %THIS IS FIGURE 4 (plotIndex 18)
    %got rid of noevolve, changed to best fastdst.
    plotIndex = 18;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','bestTreat_aveCatIV_empUpta_fastDST_',32};
    totPlotName{plotIndex} = 'summaryIntxn';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 4;
    positionMat{plotIndex} = [1:1:32];
    sumPlot = 1;
    sumYrChange{plotIndex} = [4:4:32];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %19 and 20 are for halfDOTSplus ramp up, 21 is for longDotsPlus ramp up
    plotIndex = 19;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_halfDotsPlus_',12};
    totPlotName{plotIndex} = 'summaryHalfDotsPlus';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:12];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:12];
    desiredYrs{plotIndex} = [2017, 2027];
    
    plotIndex = 20;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_',12};
    totPlotName{plotIndex} = 'summaryFullDotsPlus';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:12];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:12];
    desiredYrs{plotIndex} = [2017, 2027];
    
    plotIndex = 21;  %long dots ramp up
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','bestTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_fastDST_longDotsPlus_',12};
    totPlotName{plotIndex} = 'summaryLongDotsPlus';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:12];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:12];
    desiredYrs{plotIndex} = [2017, 2027];
    
    %Inital MDR Seed sensitivity
    plotIndex = 22;  %doubleMDRseed
    plotNameArray{plotIndex}={'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_doubleMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_doubleMDRseed_',24};
    totPlotName{plotIndex} = 'summaryDoubleMDRseed';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;  %number of years x 2 and then the number of scenarios
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 23;  %halfMDRseed
    plotNameArray{plotIndex}={'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_','bestTreat_aveCatIV_empUpta_slowDST_halfMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_halfMDRseed_',24};
    totPlotName{plotIndex} = 'summaryHalfMDRseed';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;  %number of years x 2 and then the number of scenarios
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 24;  %lessLatMDRseed
    plotNameArray{plotIndex}={'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_lessLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_lessLatMDRseed_',24};
    totPlotName{plotIndex} = 'summaryLessLatMDRseed';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;  %number of years x 2 and then the number of scenarios
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 25;  %moreLatMDRseed
    plotNameArray{plotIndex}={'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_','inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_','bestTreat_aveCatIV_empUpta_slowDST_moreLatMDRseed_','aveTreat_aveCatIV_empUpta_fastDST_moreLatMDRseed_',24};
    totPlotName{plotIndex} = 'summaryMoreLatMDRseed';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;  %number of years x 2 and then the number of scenarios
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    %THIS IS FIGURE 4 WITHOUT INTERACTION
    plotIndex = 26;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_fastDST_',24};
    totPlotName{plotIndex} = 'summaryNoIntxn';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 27;  %mdr fitness
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_','aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_',24};
    totPlotName{plotIndex} = 'summaryMDRfitness0p7';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 28;  %appendix: sensitivity analysis medium DST, no interaction
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_','aveTreat_aveCatIV_empUpta_slowDST_','bestTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_mediumDST_',24};
    totPlotName{plotIndex} = 'mediumDST';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 29;  %appendix: sensitivity analysis privateCure, no interaction
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_','aveTreat_aveCatIV_empUpta_slowDST_privateCure_','bestTreat_aveCatIV_empUpta_slowDST_privateCure_','aveTreat_aveCatIV_empUpta_fastDST_privateCure_',24};
    totPlotName{plotIndex} = 'privateCure';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 30;  %appendix: sensitivity analysis privateMDRflux0p03, no interaction
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_',24};
    totPlotName{plotIndex} = 'privateMDRflux0p03';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 31;  %appendix: sensitivity analysis privateMDRflux1p7, no interaction
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_','aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_',24};
    totPlotName{plotIndex} = 'privateMDRflux1p7';
    sizeMat{plotIndex,1} = 8; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:24];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:24];
    desiredYrs{plotIndex} = [1997, 2007, 2017, 2027];
    
    plotIndex = 32;
    plotNameArray{plotIndex}={'aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_','aveTreat_aveCatIV_empUpta_slowDST_','aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_',12};
    totPlotName{plotIndex} = 'specialDotsPlusAnalysis';
    sizeMat{plotIndex,1} = 4; sizeMat{plotIndex,2}= 3;
    positionMat{plotIndex} = [1:1:12];
    sumPlot = 1;
    sumYrChange{plotIndex} = [3:3:12];
    desiredYrs{plotIndex} = [2015,2020];
    
    % for plotNum = 11:size(plotNameArray,2)   %DEBUGGING
    for plotNum = 26   %DEBUGGING
        %for plotNum = 18  %DEBUGGING.  Run with 18 only for figure 4 of paper.  Run with 19 and 20 only for halfDOTSplus ramp up figure for SI
        aveTableTrans = zeros(sizeMat{plotNum,2}, sizeMat{plotNum,1});
        mdrtypeTypes = 2;
        if sumPlot == 1
            mdrtypeTypes = 1;
        end
        for inciPrev = 1:2
            for MDRtype = 1:mdrtypeTypes
                
                for desiredYr = desiredYrs{plotNum}
                    tempName = plotNameArray{plotNum};
                    totSubPlots = tempName{end};
                    %for each subplot, we have to grab and plot it
                    for subplotNum = 1:totSubPlots
                        desiredString = tempName{subplotNum};
                        found = 0;
                        if sumPlot == 1
                            %for making the sumGraph only
                            MDRtype = 1 ;
                            if subplotNum > (totSubPlots/2)
                                MDRtype = 2 ;
                            end
                            
                            desiredYr = desiredYrs{plotNum}(1,1) ;
                            for i = 2:size(desiredYrs{plotNum},2)
                                for j = i:size(desiredYrs{plotNum},2):size(sumYrChange{plotNum},2)
                                    if sumYrChange{plotNum}(j-1) < subplotNum && subplotNum <= sumYrChange{plotNum}(j)
                                        desiredYr = desiredYrs{plotNum}(1,i);
                                    end
                                end
                            end
                            %sprintf('subplot num %i, got desiredYr %i and MDRtype is %i',subplotNum, desiredYr,MDRtype)
                        end
                        
                        %loop over all the folders to find the right plot
                        for j = 1:size(graphFolder,2)
                            %grab the right subplot
                            if strcmp(caseStr{j},desiredString) && discontinuity{j} == desiredYr && found == 0
                                found = 1;
                                subplot(sizeMat{plotNum,1},sizeMat{plotNum,2},positionMat{plotNum}(subplotNum));
                                
                                folderNa = graphFolder{j};
                                folderNameStr = folderNa{1};
                                cd(folderNameStr);
                                if MDRtype == 1 && inciPrev == 1
                                    graphData = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
                                    graphData = [graphData(:,1)+graphData(:,4),graphData(:,2),graphData(:,3)];
                                    endTime = 516;
                                    breakTime = ((discontinuity{j}-1996)*12)+12;
                                    rowInyr = 12;
                                    mdrStr = 'MDRincidence';
                                    xvals = [1996:1/12:2038+11/12] ;
                                    ylimitMat = [0 0.00001];
                                elseif MDRtype == 2 && inciPrev == 1
                                    nonMDRinci = dlmread('TBincidence.csv', ',' ,1,0);
                                    graphData = [nonMDRinci(130:172,1), zeros(43,1), nonMDRinci(130:172,3)];
                                    endTime = size(graphData,1);
                                    breakTime = (discontinuity{j}-1996);
                                    rowInyr = 1;
                                    mdrStr = 'nonMDRincidence';
                                    xvals = [1996:1:2038];
                                    ylimitMat = [0 0.0025];
                                elseif MDRtype == 1 && inciPrev == 2
                                    xvals = [1996:1/12:2038+11/12];
                                    prevData2_tots = dlmread('monthlyActOutcomes.csv', ',' ,1,0);  %prevData2_ stuff counts transmissible MDR cases only
                                    toPlot = sum(prevData2_tots(1:516,8:10),2)./prevData2_tots(1:516,12);  %excluding col 11 since that's catIV patients
                                    graphData = [sum(prevData2_tots(1:516,8:10),2), zeros(size(toPlot,1),1),prevData2_tots(1:516,12)];
                                    breakTime = ((discontinuity{j}-1996)*12)+12;
                                    rowInyr = 12;
                                    mdrStr = 'MDRprevalence';
                                    ylimitMat = [0 0.0003];
                                elseif MDRtype == 2 && inciPrev == 2
                                    xvals = [1996:1/12:2038+11/12];
                                    prevData2_tots = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
                                    toPlot = sum(prevData2_tots(1:516,4:4),2)./prevData2_tots(1:516,12); %excluding col 5,6 since in treatment
                                    graphData = [sum(prevData2_tots(1:516,4:4),2), zeros(size(toPlot,1),1),prevData2_tots(1:516,12)];
                                    breakTime = ((discontinuity{j}-1996)*12)+12;
                                    rowInyr = 12;
                                    mdrStr = 'nonMDRprevalence';
                                    ylimitMat = [0 0.006];
                                end
                                
                                if inciPrev == 1
                                    %smooth properly
                                    %                                     timeP = {1, breakTime, breakTime, breakTime, breakTime+1, endTime};
                                    %                                     col1first = smooth((graphData(timeP{1}:timeP{2},1)+graphData(timeP{1}:timeP{2},2))./graphData(timeP{1}:timeP{2},3), smoothFactor, 'moving');
                                    %                                     col1third = smooth((graphData(timeP{5}:timeP{6},1)+graphData(timeP{5}:timeP{6},2))./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                                    %                                     col1 = [col1first;col1third];
                                    %                                     col2first = smooth(graphData(1:timeP{2},2)./graphData(1:timeP{2},3), smoothFactor, 'moving');
                                    %                                     col2third = smooth(graphData(timeP{5}:timeP{6},2)./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                                    %                                     col2 = [col2first; col2third];
                                    %                                     col3first = smooth(graphData(1:timeP{2},1)./graphData(1:timeP{2},3), smoothFactor, 'moving');
                                    %                                     col3third = smooth(graphData(timeP{5}:timeP{6},1)./graphData(timeP{5}:timeP{6},3), smoothFactor, 'moving');
                                    %                                     col3 = [col3first;col3third];
                                    col1 = smooth((graphData(1:endTime,1)+graphData(1:endTime,2))./graphData(1:endTime,3), smoothFactor, 'moving');
                                    col2 = smooth(graphData(1:endTime,2)./graphData(1:endTime,3), smoothFactor, 'moving');
                                    col3 = smooth(graphData(1:endTime,1)./graphData(1:endTime,3), smoothFactor, 'moving');
                                    toPlot = [col1, col2, col3];
                                end
                                
                                if sumPlot == 1
                                    tillStart = 1;
                                    tillEnd = 10; %5;
                                    %sprintf('start ave %i, end time %i',xvals(breakTime+floor(tillStart*rowInyr)), xvals(breakTime+floor(tillEnd*rowInyr)))
                                    %grab the average and stick it in a matrix.  start 1 year after policy starts, end 11 years after policy starts
                                    aveData = (graphData(breakTime+floor(tillStart*rowInyr):breakTime+floor(tillEnd*rowInyr),1)+graphData(breakTime+floor(tillStart*rowInyr):breakTime+floor(tillEnd*rowInyr),2))./graphData(breakTime+floor(tillStart*rowInyr):breakTime+floor(tillEnd*rowInyr),3);
                                    aveTableTrans(positionMat{plotNum}(subplotNum)) = mean(aveData,1);
                                end
                                
                                %plot the subplot
                                plot(xvals, toPlot);
                                subplotNameStr = strrep(caseStr{j}, '_', ' ');
                                underscorePositions = strfind(caseStr{j}, '_');
                                subplotNameStr1 = subplotNameStr(1:underscorePositions(3)-1);
                                subplotNameStr2 = subplotNameStr(underscorePositions(3)+1:end);
                                title({subplotNameStr1, subplotNameStr2});
                                if MDRtype == 2
                                    xlabel(num2str(toPlot(1,end)));
                                end
                                ylim(ylimitMat);
                                xlim([min(xvals) max(xvals)]);
                            end
                            if found == 1
                                break;
                            end
                        end
                    end
                    if sumPlot ~= 1
                        %save the graph
                        %   cd(plotFixerFold);
                        %   plotfixer;
                        set(gcf,'PaperPosition',[0 0  3*sizeMat{plotNum,2}  2*sizeMat{plotNum,1}])
                        typeStr = strcat(totPlotName{plotNum});
                        currentFold = pwd;
                        subOutFold = strcat(outputFolder,'\figure4\',mdrStr,'\subplots\', startStr{j});
                        print(printFormat,plotResolution,[subOutFold '/' typeStr]);
                        cd(currentFold);
                        close all
                    end
                end
                if sumPlot == 1
                    %save the graph
                    %   cd(plotFixerFold);
                    %   plotfixer;
                    set(gcf,'PaperPosition',[0 0  3*sizeMat{plotNum,2}  2*sizeMat{plotNum,1}])
                    typeStr = strcat(totPlotName{plotNum});
                    currentFold = pwd;
                    subOutFold = strcat(outputFolder,'\figure4\',mdrStr,'\subplots\', startStr{j});
                    print(printFormat,plotResolution,[subOutFold '/' typeStr]);
                    cd(currentFold);
                    close all
                    %print the average matrix
                    cd(masterFolder);
                    tableHeader = 'years on rows columns are scenarios';
                    aveTableName = strcat('aveTable', totPlotName{plotNum});
                    tablePrinter(tableHeader, aveTableTrans', aveTableName, subOutFold);
                end
            end
        end
        
        %make the comparison bar graphs
        if sumPlot == 1
            incPrevArray = {'nonMDRincidence', 'nonMDRprevalence'};
            
            for incPrev = 1:2
                dataFold = strcat('paperFigures\figure4\',incPrevArray{incPrev},'\subplots\',num2str(desiredYrs{plotNum}(1,end)));
                aveTableFolder = strcat(masterFolder,dataFold);
                cd(aveTableFolder)
                tableName = strcat(aveTableName, '.csv');
                aveTable = dlmread(tableName, ',' ,1,0);
                barPlotData{1} = aveTable(1:sizeMat{plotNum,1}/2,1:end);  %MDR
                barPlotData{2} = aveTable((sizeMat{plotNum,1}/2)+1:sizeMat{plotNum,1},1:end);  %DS
                
                %levels graph
                subplot(2,1,1)  %MDR
                bar(barPlotData{1}');
                if incPrev == 1
                    title({'MDR incidence','  (base, bestTreat, noEvolve, fastDST, accDST, bestFast, bestAcc, noEvolveFast, noEvolveAcc)'});
                elseif incPrev == 2
                    title({'MDR prevalence','  (base, bestTreat, noEvolve, fastDST, accDST, bestFast, bestAcc, noEvolveFast, noEvolveAcc)'});
                end
                subplot(2,1,2)  %DS
                bar(barPlotData{2}');
                if incPrev == 1
                    title('non-MDR incidence');
                elseif incPrev == 2
                    title('non-MDR prevalence');
                end
                legend('policy initalized 1997','policy initalized 2007', 'policy initalized 2017','policy initalized 2027', 'Location', 'SouthWestOutside');
                typeStr = strcat('aveByScenario', totPlotName{plotNum});
                print(printFormat,plotResolution,[aveTableFolder '/' typeStr]);
                close all
                
                %difference graph
                for j = 1:2
                    diffInIncidence{j} = repmat(barPlotData{j}(:,1),1,size(barPlotData{j},2)-1)-barPlotData{j}(:,2:end);
                end
                subplot(2,1,1)  %MDR
                bar(diffInIncidence{1}');
                if incPrev == 1
                    title({'Reduction in MDR incidence', '(bestTreat, noEvolve, fastDST, accDST, best fast, best acc, noEvolve fast, noEvolve acc)'})
                elseif incPrev == 2
                    title({'Reduction in MDR prevalence','(bestTreat, noEvolve, fastDST, accDST, best fast, best acc, noEvolve fast, noEvolve acc)'})
                end
                subplot(2,1,2)  %DS
                bar(diffInIncidence{2}');
                if incPrev == 1
                    title('Reduction in non-MDR incidence');
                    ylim([0 0.0003]);
                elseif incPrev == 2
                    title('Reduction in non-MDR prevalence');
                end
                legend('policy initalized 1997','policy initalized 2007','policy initalized 2017','policy initalized 2027', 'Location', 'SouthWestOutside');
                typeStr = strcat('diffInAveRateFromBase', totPlotName{plotNum});
                print(printFormat,plotResolution,[aveTableFolder '/' typeStr]);
                close all
                
                %percentage difference graph
                for j = 1:2
                    PercDiffInIncidence{j} = 100*((repmat(barPlotData{j}(:,1),1,size(barPlotData{j},2)-1) - barPlotData{j}(:,2:end)) ./  (repmat(barPlotData{j}(:,1),1,size(barPlotData{j},2)-1)));
                end
                subplot(2,1,1)  %MDR
                bar(PercDiffInIncidence{1}');
                if incPrev == 1
                    %    title({'Percentage Reduction in MDR incidence', '(bestTreat, noEvolve, fastDST, accDST, best fast, best acc, noEvolve fast, noEvolve acc)'})
                elseif incPrev == 2
                    %    title({'Percentage Reduction in MDR prevalence','(bestTreat, noEvolve, fastDST, accDST, best fast, best acc, noEvolve fast, noEvolve acc)'})
                end
                ylim([0 100]);
                set(gca,'XTick',[]);
                
                subplot(2,1,2)  %DS
                bar(PercDiffInIncidence{2}');
                if incPrev == 1
                    %    title('Percentage Reduction in non-MDR incidence');
                elseif incPrev == 2
                    %   title('Percentage Reduction in non-MDR prevalence');
                end
                set(gca,'XTick',[]);
                %legend('policy initalized 1997','policy initalized 2007','policy initalized 2017','policy initalized 2027', 'Location', 'SouthWestOutside');
                typeStr = strcat('PercDiffInAveInciRateFromBase', totPlotName{plotNum});
                print(printFormat,plotResolution,[aveTableFolder '/' typeStr]);
                close all
                
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%this does not work yet
if giantPrevGraph == 1
    %make prevalence graph
    years = [1996:5:2041];
    mdrStr = '';
    for j = 1:size(graphFolder,2)
        folderNa = graphFolder{j};
        folderNameStr = folderNa{1}(1:end);
        cd(folderNameStr);
        years = [1996:5:2041];
        prevData = dlmread('HealthOutcomes.csv', ',' ,1,0);
        mdrPrevData{j} = prevData(26:36:10);
        
        plot(years,mdrPrevData{j});
        hold on;
    end
    currentFold = pwd;
    subOutFold = strcat(outputFolder,'\figure4\',mdrStr,'\subplots\', startStr{j});
    print('-dpdf',plotResolution,[subOutFold '/' typeStr]);
    cd(currentFold);
    close all
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ratioOfTransToTx == 1
    plotResolution = '-r100';
    folderNameStr = strcat(masterFolder, 'aveTreat_aveCatIV_empUpta_slowDST_1997\2012-07-20_12-34-54');
    
    cd(folderNameStr);
    graphData = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
    endTime = 516;
    xvals = [1996:1/12:2038+11/12] ;
    
    pretoPlot = graphData(1:endTime,2)./(graphData(1:endTime,1)+graphData(1:endTime,2)) ;
    pretoPlot(pretoPlot==Inf) = max(find(pretoPlot==Inf))+1;
    toPlot = smooth(pretoPlot, 240, 'moving');
    plot(xvals(142:end-24), toPlot(142:end-24,:), 'Color', [0.7 0.5 0 ]);
    ylim([0 1]);
    ylabel({'Percentage of Incident MDR Cases','That were Transmission-Generated'});
    xlabel('Year');
    currentFold = pwd;
    cd(plotFixerFold);
    plotfixer;
    typeStr = 'PercentageTrans';
    print('-dpdf',plotResolution,[outputFolder '/' typeStr]);
    cd(currentFold);
    close all
    
    pretoPlot = graphData(1:endTime,2)./(graphData(1:endTime,1)+graphData(1:endTime,2)) ;
    pretoPlot(pretoPlot==Inf) = max(find(pretoPlot==Inf))+1;
    toPlot = smooth(pretoPlot, 240, 'moving');
    plot(xvals(142:end-24), toPlot(142:end-24,:), 'Color', [0.7 0.5 0 ]);
    hold on;
    pretoPlot = graphData(1:endTime,1)./(graphData(1:endTime,1)+graphData(1:endTime,2)) ;
    pretoPlot(pretoPlot==Inf) = max(find(pretoPlot==Inf))+1;
    toPlot = smooth(pretoPlot, 240, 'moving');
    plot(xvals(142:end-24), toPlot(142:end-24,:), 'k');
    ylim([0 1]);
    ylabel({'Percentage of Incident MDR Cases'});
    legend('Transmission-gen','tx-gen');
    xlabel('Year');
    currentFold = pwd;
    cd(plotFixerFold);
    plotfixer;
    typeStr = 'PercentageTrans';
    print('-dpdf',plotResolution,[outputFolder '/' typeStr]);
    cd(currentFold);
    close all
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if findTotInciCases == 1
    %%FINDING THE TOTAL NUMBER OF INCIDENT CASES BETWEEN 2013 AND 2038 FOR PAPER
    WHO_popProjection = [...
        2010    1224614000
        2015    1308221000
        2020    1386909000
        2025    1458958000
        2030    1523482000
        2035    1579802000
        2040    1627029000
        ];
    whoYrlyPop = interp1(WHO_popProjection(:,1), WHO_popProjection(:,2),[2013:1/12:2038+11/12]);
    cd(baseFolder);
    mdrInci = dlmread('incidenceMatrix_forMakingWHOcomparison.csv', ',' ,1,0);
    fracActivated_monthly = mdrInci(:,1)./mdrInci(:,3);
    years = [1866+1/12:1/12:2046];
    %row num 1764 is 2013 and row num 2075 is yr 2038
    totalIncidence = sum(fracActivated_monthly(1764:2075).*whoYrlyPop');
    sprintf('the total number of incident cases between 2013 and 2038 is %i', totalIncidence)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sensitivity analysis with mdr seed
if sensi_MDRseed == 1
    %get folders
    halfMDRfold = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_halfMDRSeed_1997\2012-09-20_16-37-56');
    doubleMDRfold = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_doubleMDRseed_1997\2012-09-20_14-31-42');
    lessLatfold = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_lessLatMDRseed_1997\2012-09-25_17-52-56');
    moreLatfold = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_moreLatMDRseed_1997\2012-09-25_19-55-18');
    halfDotsPlus = strcat(masterFolder, 'aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_2017\2012-08-23_17-31-34');
    longDotsPlus = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_longDotsPlus_2017\2012-09-26_21-17-37');
    
    typeName = {'base', 'half', 'double', 'lessLat', 'moreLat','halfDotsPlus','longDotsPlus'};
    typeFold = {baseFolder, halfMDRfold, doubleMDRfold, lessLatfold, moreLatfold, halfDotsPlus, longDotsPlus};
    
    %grab data
    smoothFactor = 120;
    
    for i = 1: size(typeName,2)
        cd(typeFold{i});
        fullData{i} = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        plotdata{i} = fullData{i}(:,2)./sum(fullData{i}(:,1:2),2);
        smoothed{i} = smooth(plotdata{i}, smoothFactor, 'moving');
        prev{i} = dlmread('healthOutcomes.csv', ',' ,1,0);
        prevPlot{i} = prev{i}(20:end,10);
        act{i} = dlmread('actTBoutcomes.csv', ',' ,1,0);
        actPlot{i} = sum(act{i}(:,13:16),2);
        cat4Plot{i} = act{i}(7:16,8)./sum(act{i}(7:16,5:8),2);
        cat4_smoothed{i} = cat4Plot{i};
        %cat4_smoothed{i} = smooth(cat4Plot{i}, 4, 'moving');
    end
    
    %plot
    years = [1996:1/12:2046];
    plot(years(207:517), smoothed{3}(207:517,1), 'r');
    hold on;
    plot(years(207:517), smoothed{4}(207:517,1),'m');
    hold on;
    plot(years(207:517), smoothed{1}(207:517,1),'b');
    hold on;
    plot(years(207:517), smoothed{5}(207:517,1),'g');
    hold on;
    plot(years(207:517), smoothed{2}(207:517,1),'k');
    
    %legend('Double Lat., Double Act. ','Half Lat., Double Act.','Base Case','Double Lat., Half Act.','Half Lat., Half Act.');
    legend('Lat. Prev. ~ 500, Act. Prev. ~ 6',...
        'Lat. Prev. ~ 125, Act. Prev. ~ 6',...
        'Lat. Prev. ~ 250, Act. Prev. ~ 3',...
        'Lat. Prev. ~ 500, Act. Prev. ~ 1.5',...
        'Lat. Prev. ~ 125, Act. Prev. ~ 1.5');
    ylim([0 1]);
    ylabel({'Percentage of New MDR Cases That','Are Transmission-Generated'});
    xlabel('Years');
    title({'New Transmission-generated MDR Cases by Initial MDR','Prevalence, as a Proportion of Total New MDR Cases','Legend provides 1996 MDR Prevalence as Cases out of 10000'});
    
    %print
    typeStr = 'Sensitivity_initMDRseed_percTrans';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
    
    %plot
    years = [1966:5:2046];
    plot(years(7:16), prevPlot{3}(7:16), 'r');
    hold on;
    plot(years(7:16), prevPlot{4}(7:16,1),'m');
    hold on;
    plot(years(7:16), prevPlot{1}(7:16,1),'b');
    hold on;
    plot(years(7:16), prevPlot{5}(7:16,1),'g');
    hold on;
    plot(years(7:16), prevPlot{2}(7:16,1),'k');
    
    %    legend('Double Lat., Double Act. ','Half Lat., Double Act.','Base Case','Double Lat., Half Act.','Half Lat., Half Act.','Location','SouthEast');
    legend('Lat. Prev. ~ 500, Act. Prev. ~ 6',...
        'Lat. Prev. ~ 125, Act. Prev. ~ 6',...
        'Lat. Prev. ~ 250, Act. Prev. ~ 3',...
        'Lat. Prev. ~ 500, Act. Prev. ~ 1.5',...
        'Lat. Prev. ~ 125, Act. Prev. ~ 1.5','Location','SouthEast');
    %    ylim([0 1]);
    ylabel({'Active MDR TB Prevalence'});
    xlabel('Years');
    title({'Active MDR TB prevalence, by initial MDR.','Legend provides 1996 MDR Prevalence as Cases out of 10000'});
    
    %print
    typeStr = 'Sensitivity_initMDRseed_ActPrev';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
    
    %%%%%dots plus stuff
    %plot
    years = [1996:1/12:2046];
    plot(years(207:517), smoothed{1}(207:517,1),'b');
    hold on;
    plot(years(207:517), smoothed{6}(207:517,1),'r');
    hold on;
    plot(years(207:517), smoothed{7}(207:517,1),'m');
    
    legend('Base','Half','Slow');
    ylim([0 1]);
    ylabel({'Percentage of New MDR Cases That','Are Transmission-Generated'});
    xlabel('Years');
    title({'New Transmission-generated MDR Cases by Initial MDR','Prevalence, as a Proportion of Total New MDR Cases'});
    
    %print
    typeStr = 'Sensitivity_DotsPlus_percTrans';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
    
    %plot
    years = [1966:5:2046];
    plot(years(7:16), prevPlot{1}(7:16,1),'b');
    hold on;
    plot(years(7:16), prevPlot{6}(7:16), 'r');
    hold on;
    plot(years(7:16), prevPlot{7}(7:16,1),'m');
    
    legend('Base','Half','Slow');
    %    ylim([0 1]);
    xlim([2013 2038]);
    ylabel({'Active MDR TB Prevalence as a Fraction of the Population'});
    xlabel('Years');
    title({'Active MDR TB prevalence, by DOTS-plus scenario'});
    
    %print
    typeStr = 'Sensitivity_DotsPlus_ActPrev';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
    
    %plot
    years = [1966:5:2046];
    plot(years(7:16), cat4_smoothed{1},'bo-'); %(7:16,1)
    hold on;
    plot(years(7:16), cat4_smoothed{6}, 'ro-');
    hold on;
    plot(years(7:16), cat4_smoothed{7},'mo-');
    
    legend('Base','Half','Slow','Location','NorthWest');
    %    ylim([0 1]);
    xlim([2013 2038]);
    ylabel({'Proportion of Active MDR Cases on MDR Treatment'});
    xlabel('Years');
    title({'Proportion of Active TB population in MDR treatment, by DOTS-plus Ramp Up'});
    
    %print
    typeStr = 'Sensitivity_DotsPlus_cat4';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ageStructure == 1
    cd(baseFolder);
    ageMat = dlmread('AgeStructure.csv', ',' ,1,0);
    ageMatForPlot = ageMat(:,6:end);
    
    area([1871:5:2046],ageMatForPlot);
    legend('20 and Under','21 to 40','41 to 60','61 to 80','81 and Above','Location','SouthWest');
    ylim([0 1]);
    xlabel('Year');
    ylabel('Proportion of Total Population');
    %    set(gca,'XTickLabel',{[1871:5:2046]});
    
    cd(plotFixerFold);
    plotfixer;
    typeStr = 'AgePlot';
    plotResolution = '-r300';
    print(printFormat,plotResolution,[outputFolder '/' typeStr]);
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if prevComparisonGraph == 1
    masterInputFold = masterFolder;
    
    %makes grid of graphs
    %rows: base, bestTreat, fastDST, pubGeneX
    %cols: PPM (blue base, green PPM 0.1, red PPM 1),
    %      noPrivMDR (Blue base, green 0.1, red noMDR),
    %      privGeneX (Blue base, green with private GeneX, red: privateGeneX + ppm0.1, lightBlue: privateGeneX + ppm1)
    
    inputFolderArrayListName{1} = 'base';
    inputFolderArrayList{1} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_2013';...
        };
    
    %base case
    inputFolderArrayListName{2} = 'basePPM';
    inputFolderArrayList{2} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm0p1_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        };
    inputFolderArrayListName{3} = 'baseNoPrivMDR';
    inputFolderArrayList{3} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p9_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        };
    inputFolderArrayListName{4} = 'baseGeneXpub';
    inputFolderArrayList{4} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm0p1GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        };
    
    
    %best Treat
    inputFolderArrayListName{5} = 'bestTreatPPM';
    inputFolderArrayList{5} = {'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm0p1_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        };
    inputFolderArrayListName{6} = 'bestTreatNoPrivMDR';
    inputFolderArrayList{6} = {'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_privMDR0p9_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        };
    inputFolderArrayListName{7} = 'bestTreatGeneXpub';
    inputFolderArrayList{7} = {'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm0p1GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        };
    
    %fast DST
    inputFolderArrayListName{8} = 'fastDSTppm';
    inputFolderArrayList{8} = {'aveTreat_aveCatIV_empUpta_fastDST_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_ppm0p1_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_ppm1p0_2013';...
        };
    inputFolderArrayListName{9} = 'fastDSTNoPrivMDR';
    inputFolderArrayList{9} = {'aveTreat_aveCatIV_empUpta_fastDST_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_privMDR0p9_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_privMDR0p0_2013';...
        };
    inputFolderArrayListName{10} = 'bestTGeneXpub';
    inputFolderArrayList{10} = {'aveTreat_aveCatIV_empUpta_fastDST_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_ppm0p1GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_fastDST_ppm1p0GeneXpriv_2013';...
        };
    
    
    %GeneX
    inputFolderArrayListName{11} = 'pubGeneXppm';
    inputFolderArrayList{11} = {'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm0p1_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm1p0_2013';...
        };
    inputFolderArrayListName{12} = 'pubGeneXNoPrivMDR';
    inputFolderArrayList{12} = {'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_privMDR0p9_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_privMDR0p0_2013';...
        };
    inputFolderArrayListName{13} = 'pubGeneXGeneXpub';
    inputFolderArrayList{13} = {'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm0p1GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm1p0GeneXpriv_2013';...
        };
    
    
    for graphNum = 1:4
        for z = 2:size(inputFolderArrayList,2)
            
            
            %grab the data
            years = [1870:5:2046];
            inputFolderArray = inputFolderArrayList{z};
            clearvars legStr healthTable healthTableAll_DS healthTableAll_MDR incTable incTable_DS incTable_MDR
            for i = 1:size(inputFolderArray,1)
                inputFold = strcat(masterInputFold,inputFolderArray{i});
                cd(inputFold);
                foldersStruc = dir(inputFold);
                is_dir = [foldersStruc(:).isdir]';
                foldnameArray_all = {foldersStruc(:).name}';
                foldNameArray = foldnameArray_all(is_dir);
                foldNameArray = foldNameArray(3:end); % Exclude . and ..
                cd(foldNameArray{end})
                
                healthTable{i} = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
                healthTableAll_DS(:,i) = sum(healthTable{i}(:,4:4),2)./healthTable{i}(:,12);
                healthTableAll_MDR(:,i) = sum(healthTable{i}(:,8:10),2)./healthTable{i}(:,12);
                incTable{i} = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
                incTable_DS(:,i) = incTable{i}(:,3)./incTable{i}(:,5);
                incTable_MDR(:,i) = incTable{i}(:,4)./incTable{i}(:,5);
                
                %fix legend
                legStr{i} = strrep(inputFolderArray{i}, '_', ' ');
            end
            
            %plot the data
            if graphNum == 1
                subplot(4,3,z-1);
                plot([1996:1/12:2046], healthTableAll_DS);
                typeStr = 'DS_prev';
                yRange = [0 0.0035];
                ylabel('DS Prevalence');
                xlabel('Year');
                xlim([1990 2050]);
                ylim(yRange);
                %            legend(legStr, 'Location','Best');%'NorthEast'
            elseif graphNum == 2
                subplot(4,3,z-1);
                plot([1996:1/12:2046], healthTableAll_MDR);
                typeStr = 'MDR_prev';
                yRange = [0 0.00045];
                ylabel('MDR Prevalence');
                xlabel('Year');
                xlim([1990 2050]);
                ylim(yRange);
                %            legend(legStr, 'Location','Best');%'NorthEast'
            elseif graphNum == 3
                subplot(4,3,z-1);
                plot([1996:1/12:2046], incTable_DS);
                typeStr = 'DS_inci';
                yRange = [0 0.0002];
                ylabel('DS Incidence');
                xlabel('Year');
                xlim([1990 2050]);
                ylim(yRange);
                %              legend(legStr, 'Location','Best');%'NorthEast'
            elseif graphNum == 4
                subplot(4,3,z-1);
                plot([1996:1/12:2046], incTable_MDR);
                typeStr = 'MDR_inci';
                yRange = [0 0.00002];
                ylabel('MDR Incidence');
                xlabel('Year');
                xlim([1990 2050]);
                ylim(yRange);
                %            legend(legStr, 'Location','Best');%'NorthEast'
            end
            
        end  %close z
        
        plotResolution = '-r250';
        %name = strcat(typeStr,'_',inputFolderArrayListName{z});
        print(printFormat,plotResolution,[outputFolder '/' typeStr]);
        close all;
        
    end  %close graph type
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SMDMpubPrivateGraphs == 1
    masterInputFold = masterFolder;
    
    %PUBLIC Interventions
    %inputColorList1 = [0 0 1; 1 0 1; 0 1 1];
    inputColorList1 = [0 0 1];
    inputFolderArrayListName{2} = 'base1';
    inputFolderArrayList{2} = {'aveTreat_aveCatIV_empUpta_slowDST_2013'
        };
    inputColorList{2} = inputColorList1;
    lineStyleList{2} = {'-'};
    
    i = 4;
    inputFolderArrayListName{i} = 'base2';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        };
    inputColorList{i} = inputColorList1;
    lineStyleList{i} = {'-','*','-.'};
    
    i = 5;
    inputFolderArrayListName{i} = 'base3';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        };
    inputColorList{i} = inputColorList1;
    lineStyleList{i} = {'-','*','-.'};
    
    i = 6;
    inputFolderArrayListName{i} = 'base4';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        };
    inputColorList{i} = inputColorList1;
    lineStyleList{i} = {'-','*','-.'};
    
    
    %NOT USING
    %'aveTreat_aveCatIV_empUpta_fastDST_2013'
    %'aveTreat_aveCatIV_empUpta_slowDST_GeneXpriv_2013';...
    
    
    %PRIVATE Interventions
    inputLineTypeList2 = {'-'};
    inputColorList2 = [0 0 1; 1 0.55 0; 1 0 0 ; 0 0.8 0 ; 0 0.5 0 ; 0.8 0.8 0.8 ; 0 0 0];
    
    i = 7;
    inputFolderArrayListName{i} = 'baseAllPriv1';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        };
    inputColorList{i} = [0 0 1; 1 0 0 ;   0 0.5 0  ;  0 0 0];
    legendStr{i} = {'Base','PPM 10%','PPM','MDR Red 10%', 'No Priv MDR','Priv GeneX PPM 10%','Priv GeneX PPM'};
    lineStyleList{i} = inputLineTypeList2;
    
    i = 8;
    inputFolderArrayListName{i} = 'baseAllPriv2';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        };
    inputColorList{i} = [0 0 1; 1 0 0 ;   0 0.5 0  ;  0 0 0];
    legendStr{i} = {'Base','PPM 10%','PPM','MDR Red 10%', 'No Priv MDR','Priv GeneX PPM 10%','Priv GeneX PPM'};
    lineStyleList{i} = inputLineTypeList2;
    
    i = 9;
    inputFolderArrayListName{i} = 'baseAllPriv3';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        };
    inputColorList{i} = [0 0 1; 1 0 0 ;   0 0.5 0  ;  0 0 0];
    legendStr{i} = {'Base','PPM 10%','PPM','MDR Red 10%', 'No Priv MDR','Priv GeneX PPM 10%','Priv GeneX PPM'};
    lineStyleList{i} = inputLineTypeList2;
    
    i = 10;
    inputFolderArrayListName{i} = 'baseAllPriv4';
    inputFolderArrayList{i} = {'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        };
    inputColorList{i} = [0 0 1; 1 0 0 ;   0 0.5 0  ;  0 0 0];
    legendStr{i} = {'Base','PPM 10%','PPM','MDR Red 10%', 'No Priv MDR','Priv GeneX PPM 10%','Priv GeneX PPM'};
    lineStyleList{i} = inputLineTypeList2;
    
    %Synergies
    inputColorList3 = [0 0 1; 1 0 0 ;   0 0.5 0  ;  0 0 0];
    lineStyleList3 = {'-','*','-.'};
    
    i = 11;
    inputFolderArrayListName{i} = 'synergies1';
    inputFolderArrayList{i} = {...
        'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        };
    inputColorList{i} = inputColorList3;
    lineStyleList{i} = lineStyleList3;
    
    i = 12;
    inputFolderArrayListName{i} = 'synergies2';
    inputFolderArrayList{i} = {...
        'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        };
    inputColorList{i} = inputColorList3;
    lineStyleList{i} = lineStyleList3;
    
    i = 13;
    inputFolderArrayListName{i} = 'synergies3';
    inputFolderArrayList{i} = {...
        'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        };
    inputColorList{i} = inputColorList3;
    lineStyleList{i} = lineStyleList3;
    
    i = 14;
    inputFolderArrayListName{i} = 'synergies4';
    inputFolderArrayList{i} = {...
        'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        };
    inputColorList{i} = inputColorList3;
    lineStyleList{i} = lineStyleList3;
    
    i = 15;
    inputFolderArrayListName{i} = 'synergies5';
    inputFolderArrayList{i} = {...
        'aveTreat_aveCatIV_empUpta_slowDST_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_privMDR0p0_2013';...
        'bestTreat_aveCatIV_empUpta_slowDST_ppm1p0GeneXpriv_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm1p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_privMDR0p0_2013';...
        'aveTreat_aveCatIV_empUpta_slowDST_GeneXpub_ppm1p0GeneXpriv_2013';...
        };
    inputColorList{i} = inputColorList3;
    lineStyleList{i} = lineStyleList3;
    
    
    %for graphNum = 1:4
    for graphNum = [2]
        %  for z = 2:size(inputFolderArrayList,2)
        for z = [4:1:15]
            
            %grab the data
            years = [1870:5:2046];
            inputFolderArray = inputFolderArrayList{z};
            clearvars legStr healthTable healthTableAll_DS healthTableAll_MDR incTable incTable_DS incTable_MDR smoothed_healthTableAll_DS smoothed_healthTableAll_MDR smoothed_incTable_DS smoothed_incTable_MDR
            for i = 1:size(inputFolderArray,1)
                inputFold = strcat(masterInputFold,inputFolderArray{i});
                cd(inputFold);
                foldersStruc = dir(inputFold);
                is_dir = [foldersStruc(:).isdir]';
                foldnameArray_all = {foldersStruc(:).name}';
                foldNameArray = foldnameArray_all(is_dir);
                foldNameArray = foldNameArray(3:end); % Exclude . and ..
                cd(foldNameArray{end})
                
                healthTable{i} = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
                healthTableAll_DS(:,i) = sum(healthTable{i}(:,4:4),2)./healthTable{i}(:,12);
                healthTableAll_MDR(:,i) = sum(healthTable{i}(:,8:10),2)./healthTable{i}(:,12);
                incTable{i} = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
                incTable_DS(:,i) = incTable{i}(:,3)./incTable{i}(:,5);
                incTable_MDR(:,i) = incTable{i}(:,4)./incTable{i}(:,5);
                
                %smooth
                smoothed_healthTableAll_DS(:,i) = smooth((100000*healthTableAll_DS(:,i)),12,'moving');
                smoothed_healthTableAll_MDR(:,i) = smooth((100000*healthTableAll_MDR(:,i)),36,'moving');
                smoothed_incTable_DS(:,i) = smooth((100000*incTable_DS(:,i)),12,'moving');
                smoothed_incTable_MDR(:,i) = smooth((100000*incTable_MDR(:,i)),36,'moving');
                
                %fix legend
                legStr{i} = strrep(inputFolderArray{i}, '_', ' ');
            end
            
            %plot the data
            timeVec = [1996:1/12:2046];
            %xlimVec = [1996 2050];
            %xlimVec = [1996 2023];
            xlimVec = [2011 2023];
            %             co = get(0,'defaultaxescolororder');  %Defaults
            %             lso = get(0,'defaultaxeslinestyleorder');
            %             set(0,'defaultaxescolororder',co); %Restore Defaults
            %             set(0,'defaultaxeslinestyleorder',lso);
            
            set(gca, 'ColorOrder', inputColorList{z}, 'LineStyleOrder', lineStyleList{z}, 'NextPlot', 'replacechildren');   % Change to new colors.
            %             co = get(gca,'ColorOrder') % Verify it changed
            %             lsa = get(gca,'LineStyleOrder') % Verify it changed
            
            %co = get(gca,'DefaultAxesColorOrder')  %change color order back to default
            
            if graphNum == 1
                %subplot(4,3,z-1);
                plot(timeVec(1:end-6), smoothed_healthTableAll_DS(1:end-6,:));
                typeStr = 'DS_prev';
                yRange = [0 0.0035];
                ylabel('Infectious DS Prevalence (Cases in 100,000)');
                xlabel('Year');
                xlim(xlimVec);
                ylim(yRange);
                %    legend(legStr, 'Location','Best');%'NorthEast'  'Best'
            elseif graphNum == 2
                %subplot(4,3,z-1);
                plot(timeVec(1:end-6), smoothed_healthTableAll_MDR(1:end-6,:)); %start at 121 for starting at 2006
                typeStr = 'MDR_prev';
                %                 yRange = [0 0.00045];
                yRange = [0 30];
                ylabel({'Infectious MDR Prevalence', '(Cases in 100,000)'});
                xlabel('Year');
                xlim(xlimVec);
                ylim(yRange);
                %    legend(legStr, 'Location','SouthEast');%'NorthEast''Best'
                hold on;
                plot([2012 2012], [yRange],'k--');
                set(gca,'XTick',[2010:5:2023])
            elseif graphNum == 3
                %subplot(4,3,z-1);
                plot(timeVec(1:end-6), smoothed_incTable_DS(1:end-6,:));
                typeStr = 'DS_inci';
                yRange = [0 0.0002];
                ylabel('Infectious DS Incidence (Cases in 100,000)');
                xlabel('Year');
                xlim(xlimVec);
                ylim(yRange);
                %              legend(legStr, 'Location','Best');%'NorthEast'
            elseif graphNum == 4
                %subplot(4,3,z-1);
                plot(timeVec(1:end-6), smoothed_incTable_MDR(1:end-6,:));
                typeStr = 'MDR_inci';
                yRange = [0 0.00002];
                ylabel('Infectious MDR Incidence (Cases in 100,000)');
                xlabel('Year');
                xlim(xlimVec);
                ylim(yRange);
                %            legend(legStr, 'Location','Best');%'NorthEast'
            end
            
            cd(plotFixerFold);
            plotfixerPpt;
            plotResolution = '-r250';
            name = strcat(typeStr,'_',inputFolderArrayListName{z});
            print(printFormat,plotResolution,[outputFolder '/' name]);
            close all;
        end  %close z
        
        
    end  %close graph type
    
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if prevComparisonTables ==1
    %makes a table of DS and MDR prevalence and incidence for the specified year\
    %the rows are base, best, and fast RNTCP treatment
    %the columns are the private clinic scenarios specifed
    
    %folder names are aveTreat_aveCatIV_empUpta_fastDST_ppm0p1_2013, for private scen and aveTreat_aveCatIV_empUpta_fastDST_2013 for base
    masterInputFold = masterFolder;
    cd(masterInputFold);
    cols = {'','_ppm0p1','_ppm1p0', '_privMDR0p9', '_privMDR0p0', '_GeneXpriv','_ppm0p1GeneXpriv','_ppm1p0GeneXpriv'};  %base will be first col, then the ones specified
    rows = {'aveTreat_aveCatIV_empUpta_slowDST', 'bestTreat_aveCatIV_empUpta_slowDST','aveTreat_aveCatIV_empUpta_fastDST','aveTreat_aveCatIV_empUpta_slowDST_GeneXpub'};
    yrspec = '2023';
    
    for i = 1:size(rows,2)
        for j = 1:size(cols,2)
            
            %specify which folder to use
            foldName{i,j} = strcat(rows{i},cols{j},'_2013');
            fullFoldName = strcat(masterInputFold,foldName{i,j});
            cd(fullFoldName);
            
            %grab input data
            foldersStruc = dir(fullFoldName);
            is_dir = [foldersStruc(:).isdir]';
            foldnameArray_all = {foldersStruc(:).name}';
            foldNameArray = foldnameArray_all(is_dir);
            foldNameArray = foldNameArray(3:end); % Exclude . and ..
            csvFolder = strcat(fullFoldName,'\',foldNameArray{end});
            
            %make cases csv tables
            cd(masterInputFold)
            casesAvertedMaker(masterInputFold, 0, '', strcat(foldName{i,j},'\', foldNameArray{end}),2023, 'tot');
            cd(csvFolder)
            
            prevTable = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
            incTable = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
            casesTable = dlmread(strcat('casesTotal_from1996_',yrspec,'_tot.csv'), ',' ,1,0);
            
            %make data for table
            prev_DS(i,j)  = mean(sum(prevTable(313:348,4:4),2))./mean(prevTable(313:348,12));
            prev_MDR(i,j) = mean(sum(prevTable(313:348,8:10),2))./mean(prevTable(313:348,12));
            
            inc_DS(i,j)  = 12*mean(incTable(313:348,3)./mean(incTable(313:348,5)));   %monthly incidence rate (averaged over 3 years) times 12 for annual inc 2023
            inc_MDR(i,j) = 12*mean(incTable(313:348,4)./mean(incTable(313:348,5)));
            
            cases_DS(i,j)  = casesTable(1,3);
            cases_MDR(i,j) = casesTable(1,4);
            
            if i == 1 && j == 1
                %get base case 2013 values for making comparisons later
                prev_DS(size(rows,2)+1,size(cols,2)+1)  = mean(sum(prevTable(193:228,4:4),2))./mean(prevTable(193:228,12));  %want to use year 2013, so average between rows 31 and 32 on healthOutcomes table
                prev_MDR(size(rows,2)+1,size(cols,2)+1) = mean(sum(prevTable(193:228,8:10),2))./mean(prevTable(193:228,12));
                inc_DS(size(rows,2)+1,size(cols,2)+1)  = 12*mean(incTable(193:228,3)./mean(incTable(193:228,5)));   %monthly incidence rate (averaged over 3 years) times 12 for annual inc 2013
                inc_MDR(size(rows,2)+1,size(cols,2)+1) = 12*mean(incTable(193:228,4)./mean(incTable(193:228,5)));
                
                cd(masterInputFold)
                casesAvertedMaker(masterInputFold, 0, '', strcat(foldName{i,j},'\', foldNameArray{end}),2013, 'tot');
                cd(csvFolder)
                casesTable = dlmread(strcat('casesTotal_from1996_','2013','_tot.csv'), ',' ,1,0);
                cases_DS(size(rows,2)+1,size(cols,2)+1) = casesTable(1,3);
                cases_MDR(size(rows,2)+1,size(cols,2)+1) = casesTable(1,4);
            end
        end
    end
    
    cd(masterInputFold);
    tableHeader = 'In 2023. Base in 2013 bottom right. rows: base best and fast RNTCP treatment. Col: base, ppm0p1, ppm1p0, privMDR0, privMDR0p5,GeneXpriv,ppm1p0GeneXpriv';
    tablePrinter(tableHeader, prev_DS, strcat('prev_',yrspec,'_DS'), outputFolder);
    tablePrinter(tableHeader, prev_MDR, strcat('prev_',yrspec,'_MDR'), outputFolder);
    tablePrinter(tableHeader, inc_DS, strcat('inc_',yrspec,'_DS'), outputFolder);
    tablePrinter(tableHeader, inc_MDR, strcat('inc_',yrspec,'_MDR'), outputFolder);
    tablePrinter(tableHeader, cases_DS, strcat('casesSince1996_',yrspec,'_DS'), outputFolder);
    tablePrinter(tableHeader, cases_MDR, strcat('casesSince1996_',yrspec,'_MDR'), outputFolder);
end


%=============================

if sensitivityAllMDR == 1
    
    %may need to make the graph and then manually save to get the legend to
    %not cover stuff up
    
    plotNameArray1 = {...
        'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_privateCure_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_1997';...
        %'aveTreat_aveCatIV_empUpta_slowDST_vyn_1997';...
        };
    
    plotNameArray2 = {
        'bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_longDotsPlus_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_privateCure_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_2017';...
        'bestTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_2017';...
        };
    
    plotNameArray3 = {...
        'aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux1p7_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_halfDotsPlus_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_longDotsPlus_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_privateCure_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_mdrFit0p7_2017';...
        'aveTreat_aveCatIV_empUpta_fastDST_privatMDRflux0p03_2017';...
        };
    
    plotNameArray = [plotNameArray1,plotNameArray2,plotNameArray3];
    
    for j = 1:size(plotNameArray,2)
        for i = 1:size(plotNameArray,1)
            folder = strcat(masterFolder,plotNameArray{i,j});
            cd(folder)
            foldersStruc = dir(folder);
            is_dir = [foldersStruc(:).isdir]';
            foldnameArray_all = {foldersStruc(:).name}';
            foldNameArray = foldnameArray_all(is_dir);
            foldNameArray = foldNameArray(3:end); % Exclude . and ..
            cd(foldNameArray{1})
            
            healthOutcomes = dlmread('healthOutcomes.csv', ',' ,1,0);
            actTBoutcomes = dlmread('actTBoutcomes.csv', ',' ,1,0);
            infectiousMDR = sum(actTBoutcomes(7:end,5:7),2);
            totalNumPpl = sum(healthOutcomes(26:end,1:5),2);
            
            MDRlevel{j}(:,i) = infectiousMDR./totalNumPpl;
        end
        
        NumMDRpplOutOfHunThou{j} = MDRlevel{j}*100000;
        
        if j == 1
            titleStr = 'sensitivityMDRsummary_base';
        elseif j == 2
            titleStr = 'sensitivityMDRsummary_best2017';
        elseif j == 3
            titleStr = 'sensitivityMDRsummary_fast2017';
        end
        
        %plot
        time = [ 1996        2001        2006        2011        2017        2021        2026        2031 2036];
        plot(time,NumMDRpplOutOfHunThou{j}(1:9,:));
        hold on;
        plot([2006 2006], [0 40], 'k-.');
        hold on;
        plot([2013 2013], [0 40], 'k:');
        
        if j == 1
            legend('Scenario 1','Scenario 2','Scenario 3','Base Case','Scenario 4','Scenario 5','Scenario 6','DOTS-Plus Begins','Current Time','Location','SouthEast');
            
        else
            hold on;
            plot([2017 2017], [0 40], 'k--');
            %             if j ==2  %run with this uncommented to get the legend;
            %                 legend('Scenario 1','Scenario 2','Scenario 3','Base Case','Scenario 4','Scenario 5','Scenario 6','DOTS-Plus Begins','Current Time','Policy Initiation','Location','SouthEast');
            %                 titleStr = strcat(titleStr,'wLegend');
            %             elseif j ==3
            %                 legend('Scenario 1','Scenario 2','Scenario 3','Base Case','Scenario 4','Scenario 5','Scenario 6','DOTS-Plus Begins','Current Time','Policy Initiation','Location','NorthEast');
            %                 titleStr = strcat(titleStr,'wLegend');
            %             end
        end
        % legend('1.7x Private MDR','Incomplete Dots+ Ramp Up','Long Dots+ Ramp Up','Base Case','Private Cures Non-MDR TB','0.7 MDR Fitness','0.03x Private MDR','Current Time','Location','SouthEast');
        ylabel('Number of Infectious MDR TB Cases Per 100,000 People');
        xlabel('Time');
        title(' MDR TB Prevalence Under Various Sensitivity Scenarios');
        xlim([1996,2036]);
        ylim([0 40])
        cd(plotFixerFold);
        plotFixerPNAS;
        
        plotResolution = '-r170';
        print(printFormat,plotResolution,[outputFolder '/' titleStr]);
        close all;
    end
    
end


%======================================

if make2013PropTransGenMDR == 1
    
    
    plotNameArray = {...
        'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux1p7_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_halfDotsPlus_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_longDotsPlus_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_privateCure_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_mdrFit0p7_1997';...
        'aveTreat_aveCatIV_empUpta_slowDST_privatMDRflux0p03_1997';...
        %'aveTreat_aveCatIV_empUpta_slowDST_vyn_1997';...
        };
    
    
    for j = 1:size(plotNameArray,1)
        folder = strcat(masterFolder,plotNameArray{j});
        cd(folder)
        foldersStruc = dir(folder);
        is_dir = [foldersStruc(:).isdir]';
        foldnameArray_all = {foldersStruc(:).name}';
        foldNameArray = foldnameArray_all(is_dir);
        foldNameArray = foldNameArray(3:end); % Exclude . and ..
        cd(foldNameArray{1})
        
        healthOutcomes = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
        transMDR_2013{j} = mean(healthOutcomes(193:217,2)./sum(healthOutcomes(193:217,[1,2,4]),2));
    end
    transMDR_2013
end



%======================================
if TB_MACtargets == 1
    
    %  for 2000 - 2025
    % TB notification (n/100k/year) (all ages, all types, all sectors)
    % TB prevalence (n/100k/) (all ages, all types)
    % TB mortality (n/100k/year) (all ages, all types, all sectors)
    % TB disease incidence (n/100k/year)  (all ages, all types)
    % TB disease prevalence by MDR status (n/100k/) (smear positive, adults, all treatment histories)
    % Proportion of all incident TB disease that is in <15 year olds (all types, all sectors)
    % Proportion of total population latently infected with M.tb
    IndianPop = [...
        1950    371857000   ;...
        1955    406374000   ;...
        1960    447844000   ;...
        1965    496400000   ;...
        1970    553874000   ;...
        1975    622097000   ;...
        1980    700059000   ;...
        1985    784491000   ;...
        1990    873785000   ;...
        1995    964486000   ;...
        2000    1053898000  ;...
        2005    1140043000  ;...
        2010    1224614000  ;...
        2015    1308221000  ;...
        2020    1386909000  ;...
        2025    1458958000  ;...
        2030    1523482000  ;...
        2035    1579802000  ;...
        2040    1627029000  ;...
        2045    1664519000  ;...
        2050    1692008000  ;...
        ];
    TBparams.IndianPop = interp1(IndianPop(:,1), IndianPop(:,2), [1990:1/12:2046]);
    popYear = [1990:1/12:2046]';
    
    for i = 1:6
        toPlot{i} = [];
    end
    
    %initalize
    TB_MAC_epi_results_all = []; TB_MAC_daly1_results_all = []; TB_MAC_daly2_results_all= []; TB_MAC_daly3_results_all = []; TB_MAC_econ_results_all = [];
    
    for folderNum = 1 : size(TB_MACtargetFolds,2)
        %grab folder
        cd(masterFolder);
        folderNameStartsWith = TB_MACtargetFolderStr;  %this is the first letter of the runs 'p01'
        mustHaveCSV = 'diagnosedPpl.csv';
        currentFold = findMostRecentLegitFolder(TB_MACtargetFolds{folderNum}, folderNameStartsWith, mustHaveCSV)
        %currentFold = fullfile(masterFolder, TB_MACtargetFolds{folderNum})
        cd(currentFold);
        
        %grab data
        diagnosedPpl = dlmread('diagnosedPpl.csv', ',' ,1560,0);  %includes burn in
        notification_TBmac = dlmread('notification_TBmac.csv', ',' ,1560,0);  %includes burn in
        monthlyActiveOutcomes = dlmread('monthlyActOutcomes.csv', ',' ,1,0);  %post burn in
        deathsCounterMat = dlmread('deathsCounterMat.csv', ',' ,1560,0); %includes burn in
        deathsTBAgesMat = dlmread('deathTBAgesMat_TBmac.csv', ',' ,1,0);
        deathAgesMat = dlmread('deathAgesMat_TBmac.csv', ',' ,1,0);
        latentToActive = dlmread('latentToActive.csv', ',' ,1560,0);  %includes burn in
        numMDRinc = dlmread('numMDRactivations_postBurnIn.csv',',',1,0);
        numPplInDots = dlmread('DOTSperformanceMeas.csv',',',1,0);
        latIncidence = dlmread('DSmdrIncidence_postBurnIn.csv',',',1,0);
        pastTreat = dlmread('pastTreatmentCounts.csv',',',1,0);  %active TB only, cols are pastTrt = 0, 1, 2, 3, and any past private
        private = dlmread('privateCounterMat.csv',',',1560,0); % col 4 diagnosed private this month %includes burn in
        aliveAgesMat = dlmread('aliveAgesMat_TBmac.csv', ',' ,1,0);
        aliveTBAgesMat = dlmread('aliveTBAgesMat_TBmac.csv', ',' ,1,0);   %ages of TB ppl
        diagnosisCounterMat = dlmread('diagnosedPpl.csv', ',' ,1560,0);
        privateCounterMat  = dlmread('privateCounterMat.csv', ',' ,1560,0);
        ACFnumScreened_TBmac = dlmread('numScreened_TBmac.csv', ',' ,1560,0);  %for intervention 4 and 5
        
        %grab the data we want
        numPpl = monthlyActiveOutcomes(:,12);
        
        notificationCATI =  diagnosedPpl(:,6);
        notificationCATII = diagnosedPpl(:,7);
        notificationCATIV = diagnosedPpl(:,5);
        notification =  notificationCATI +  notificationCATII;
        
        notifications_Revised = notification_TBmac(:,1)+ notification_TBmac(:,4); %Number of diagnosed in cat I IIand private not double counting in Col 1
        notifications_Revised_adultsOnly = notification_TBmac(:,5)+ notification_TBmac(:,6); %Number of diagnosed in cat I IIand private not double counting in Col 1
        
        TBmort =  sum(deathsCounterMat(:,[3,6,8,10]),2);
        TBmort_everTreatment =  sum(deathsCounterMat(:,[22,21,6]),2); % deathsCounterMat (22 ever private) (21 past RNTCP) (6 Cat I)
        TBmortExcess = TBmort - deathsCounterMat(:,23);
        
        TBmort_under15 =  sum(deathsTBAgesMat(:,[1,2]),2);
        TBmort_allAges =  sum(deathsTBAgesMat,2);  %check to see if the same as TBmort
        TBmort_over15 = TBmort_allAges - TBmort_under15;
        
        allMort_agesMat = [deathAgesMat(:,1),sum(deathAgesMat(:,[2,3]),2), deathAgesMat(:,[4:10]), sum(deathAgesMat(:,[11,12]),2)];
        
        alive_under15 = sum(aliveAgesMat(:,[1,2]),2); %all alive
        alive_all = sum(aliveAgesMat,2); %all alive
        
        actTBalive_over15 = sum(aliveTBAgesMat,2) - sum(aliveTBAgesMat(:,[1,2]),2); %all TB adults
        
        alive_ages = [aliveAgesMat(:,1),sum(aliveAgesMat(:,[2,3]),2), aliveAgesMat(:,[4:10]), sum(aliveAgesMat(:,[11,12]),2)];
        actTBalive_ages = [aliveTBAgesMat(:,1),sum(aliveTBAgesMat(:,[2,3]),2), aliveTBAgesMat(:,[4:10]), sum(aliveTBAgesMat(:,[11,12]),2)];
        noActTBalive_ages = alive_ages - actTBalive_ages;
        
        TBinc = sum(latentToActive(:,[3,5]),2);  % latentToActive col 3 and 5
        MDRinc = sum(numMDRinc(:,1),2);  %MDR incidence for evolved, transmitted, and private
        fracIncIsMDR = MDRinc;
        TBincUnder15 = latentToActive(:,7);
        fracFastAct = latentToActive(:,3)./ TBinc;
        
        numPplCAT12 = sum(numPplInDots(:,[1,5]),2);
        numPplCat4 = numPplInDots(:,9);
        
        % treatment success rates
        notificationCATI_offset =  [zeros(6,1);diagnosedPpl(1:end-6,6)]; %col 6 enter CATI 6 months ago
        notificationCATII_offset = [zeros(8,1);diagnosedPpl(1:end-8,7)]; %col 7 enter CATII 8 months ago
        notificationCATIV_offset = [zeros(24,1);diagnosedPpl(1:end-24,5)]; %col 5enter CATII 24 months ago
        trtSuccessCatI = numPplInDots(:,3) ./ notificationCATI_offset; %numPplInDots col 3 is catI graduating, col 4 is catI cured.  Divide by total people entering CATI 6 months ago
        trtSuccessCatII = numPplInDots(:,7) ./ notificationCATII_offset; %numPplInDots col 7 is catII graduating, col 8 is catII cured.
        trtSuccessDOTS_cat12 = (numPplInDots(:,3) + numPplInDots(:,7)) ./ (notificationCATI_offset + notificationCATII_offset);
        trtSuccessCatIV = numPplInDots(:,10) ./ notificationCATIV_offset; %numPplInDots col 10 is catIV cured
        
        annRiskofInfection = sum(latIncidence(:,[1,6]),2);  %n infections for 100000 population
        
        TBdsPrev = sum(  monthlyActiveOutcomes(:,4:7)  ,2); % monthlyActiveOutcomes col 4 to 7
        TBmdrPrev = sum(  monthlyActiveOutcomes(:,8:11)  ,2); %monthlyActiveOutcomes col 8 to 11
        TBprev =  TBdsPrev+ TBmdrPrev;
        MDRinCatI =  monthlyActiveOutcomes(:,9) ./ (monthlyActiveOutcomes(:,5) + monthlyActiveOutcomes(:,9)) ;
        MDRinCatII =  monthlyActiveOutcomes(:,10) ./ (monthlyActiveOutcomes(:,6) + monthlyActiveOutcomes(:,10));
        untreateedTB_personMon = sum(  monthlyActiveOutcomes(:,[4,8])  ,2);  %untreated
        catI_personMon = sum(  monthlyActiveOutcomes(:,[5,9])  ,2);  %in catI
        catII_personMon = sum(  monthlyActiveOutcomes(:,[6,10])  ,2);  %in catII
        catIV_personMon = sum(  monthlyActiveOutcomes(:,[7,11])  ,2);  %in catIV
        
        fracExposedToDOTScare = (sum(pastTreat(:,[2,3,4]),2)) ./ TBprev;  %catI is marked as prevTrted already plus previously treated
        fracExposedToPrivatecare = (private(:,4) + pastTreat(:,5)) ./ TBprev;  %in private and have exposure to private
        inPrivate = private(:, 4);
        inPrivate_privRetrt = inPrivate - notification_TBmac(:,7);  %private trt experienced
        inPrivate_trtNaive = notification_TBmac(:,7);  %private trt naive
        
        aveTime_diagRNTCP = diagnosisCounterMat(:,8);  %ave time to diagnosis in RNTCP
        aveTime_diagPrivate = privateCounterMat(:,5);  %ave time to diagnosis in private
        ratioDiagTime_RNTCPtoPriv = aveTime_diagRNTCP ./ aveTime_diagPrivate;
        RNTCPdiagWeight = notification_TBmac(:,1) ./ (notification_TBmac(:,1)+ notification_TBmac(:,4));
        aveTimeToDiag = (RNTCPdiagWeight.*aveTime_diagRNTCP) + ((1-RNTCPdiagWeight).*aveTime_diagPrivate);
        
        TBdsLatPrev = monthlyActiveOutcomes(:,2) ; %monthlyActiveOutcomes col 2
        TBmdrLatPrev =   monthlyActiveOutcomes(:,3) ; %monthlyActiveOutcomes col 3
        TBlatPrev = TBdsLatPrev +  TBmdrLatPrev;
        
        numberScreenedViaACF = ACFnumScreened_TBmac(:,1);
        diagViaACF = ACFnumScreened_TBmac(:,2);
        screenedLTBI_ACF =  ACFnumScreened_TBmac(:,3);  % 7a. Total number of persons screened for latent TB
        trtLTBI_ACF =  ACFnumScreened_TBmac(:,4);  %7b. TLTBI treatment volume
        
        time = [1996:1/12: 2036]';
        
        %inflate to real population size (just the number of people-months in dots)
        actualPop = TBparams.IndianPop(find(popYear == 1996):find(popYear == 2036));
        numPplCAT12 = actualPop'.*(numPplCAT12./numPpl);
        numPplCat4 = actualPop'.*(numPplCat4./numPpl);
        
        %collapse to yearly and output.  ONLY EVER ADD THINGS TO THE END OF THIS LIST.  WRITE NAME IN THE INDEX LIST BELOW.
        needsSumming = [notificationCATI,notificationCATII, notificationCATIV, notification, TBmort,TBmort_everTreatment,TBinc,...
            fracIncIsMDR, TBincUnder15, numPplCAT12, numPplCat4, annRiskofInfection, TBmortExcess,TBmort_under15,TBmort_allAges,...
            notifications_Revised, notifications_Revised_adultsOnly, untreateedTB_personMon, numberScreenedViaACF, diagViaACF,...
            screenedLTBI_ACF trtLTBI_ACF, catI_personMon, catII_personMon, catIV_personMon, inPrivate, inPrivate_privRetrt,  inPrivate_trtNaive];
        counter = 1;
        for yrIndex = [0:12:470]
            %index1,        2,              3,               4,               5,              6,        7,              8,       9,             10,          11,           12,           13,            14          15,               16,                17,                  18 ,                             19,                           20,                21,         22,         23,              24,            25,              26,             27,           28,                 29
            % numPpl, notificationCATI,notificationCATII, notificationCATIV, notification, TBmort,TBmort_everTreatment,TBinc,fracIncIsMDR, TBincUnder15, numPplCAT12, numPplCat4, annRiskofInfection, TBmortExcess, TBmort_under15,TBmort_allAges , notifications_Revised, notifications_Revised_adultsOnly, untreateedTB_personMon, numberScreenedViaACF, diagViaACF,screenedLTBI_ACF,trtLTBI_ACF, catI_personMon, catII_personMon, catIV_personMon, inPrivate, inPrivate_privRetrt,  inPrivate_trtNaive
            outMat1(counter, 1:size(needsSumming,2)+1) = [numPpl(yrIndex+6,1), sum(needsSumming(yrIndex+1:yrIndex+12, :) ,1) ];
            outMat2(counter, 1:14) = [TBdsPrev(yrIndex+6,1), TBmdrPrev(yrIndex+6,1),TBprev(yrIndex+6,1), TBdsLatPrev(yrIndex+6,1), TBmdrLatPrev(yrIndex+6,1),TBlatPrev(yrIndex+6,1),alive_under15(yrIndex+6,1), alive_all(yrIndex+6,1), actTBalive_over15(yrIndex+6,1), MDRinCatI(yrIndex+6,1), MDRinCatII(yrIndex+6,1), fracFastAct(yrIndex+6,1), ratioDiagTime_RNTCPtoPriv(yrIndex+6,1),aveTimeToDiag(yrIndex+6,1) ];
            outMat3(counter, 1:6) = [mean(trtSuccessCatI(yrIndex+1:yrIndex+12,1)), mean(trtSuccessCatII(yrIndex+1:yrIndex+12,1)), mean(trtSuccessDOTS_cat12(yrIndex+1:yrIndex+12,1)), mean(trtSuccessCatIV(yrIndex+1:yrIndex+12,1)), fracExposedToDOTScare(yrIndex+6,1), fracExposedToPrivatecare(yrIndex+6,1)];
            outMat4(counter, 1:30) = [allMort_agesMat(yrIndex+6,:),actTBalive_ages(yrIndex+6,:), noActTBalive_ages(yrIndex+6,:)];
            counter = counter + 1;
        end
        outMat1(:, 9) = outMat1(:, 9) ./ outMat1(:, 8);  %fracIncIsMDR is actually MDRinc / TBinc
        outMat1(:, 10) = outMat1(:, 10) ./ outMat1(:, 8);  %fraction of TBincUnder15 is TBincUnder15 / TBinc
        outMat1(:, 15) = outMat1(:, 15) ./ outMat1(:, 16);  %fraction of TBincUnder15 is TBincUnder15 / TBinc
        outMat2(:,7)  = outMat2(:,7) ./ outMat2(:,8);  %fraction alive under 15
        % outMat3  %debugging -- needed to see this to get the right intervention 1 and 2 parameters
        
        outMatFinal = [ [1996:1:2035]', ... %time
            outMat1(:, 1),...    %pop
            (outMat1(:, [2:7, 14,17, 18])./ repmat(outMat1(:, 1), 1, 9) *100000) ,... % notificationCATI,notificationCATII, notificationCATIV, notification, TBmort,TBmort_everTreatment, TBmortExcess, %notifications_revised, notifications_Revised_adultsOnly
            (outMat1(:, 8)./ repmat(outMat1(:, 1), 1, 1) *100000) ,...   %TBinc
            outMat1(:, 9:12),...                                         % fracIncIsMDR, fracTBincUnder15, numPplCAT12, numPplCat4
            outMat1(:, 13)./ outMat1(:, 1),...                           % annRiskofInfection
            (outMat2(:, 1:6)./ repmat(outMat1(:, 1), 1, 6) *100000), ... %TBdsPrev, TBmdrPrev,TBprev, TBdsLatPrev, TBmdrLatPrev,TBlatPrev
            outMat3(:, 1:6),...                                         %trtSuccessCatI, trtSuccessCatII, trtSuccessDOTS_cat12, trtSuccessCatIV, fracExposedToDOTScare, fracExposedToPrivatecare
            outMat1(:, 15),...                                          %frac TBmort_under15
            outMat2(:,7)                                               %frac_aliveUnder15
            ];
        header1 = 'Year,  pop, notificationCATI,notificationCATII, notificationCATIV, DOTSnotification, TBmort,TBmort_everTreatment,TBmortExcess, notifications_revised, notifications_Revised_adultsOnly, TBinc, fracIncIsMDR, fracTBincUnder15,numPplCAT12, numPplCat4, annRiskofInfection,';
        header2 = ' TBdsPrev, TBmdrPrev,TBprev, TBdsLatPrev, TBmdrLatPrev,TBlatPrev,trtSuccessCatI, trtSuccessCatII, trtSuccessDOTS_cat12, trtSuccessCatIV, fracExposedToDOTScare, fracExposedToPrivatecare,frac_TBmort_under15,frac_aliveUnder15';
        cd(masterFolder);
        tablePrinter(strcat(header1,header2), outMatFinal, 'TB_MAC_Output', currentFold);
        
        %make a formatted table for easy copy-pasting
        formattedOutMatFinal = [ [1996:1:2035]', ...                        % time
            outMat1(:, 1),...                                               % pop
            9912399*ones(size(outMat1(:, 1))),...							% blank
            sum(outMat1(:, 2:4)./ repmat(outMat1(:, 1), 1, 3) *100000,2) ,... % total cat I, II, and IV notifications
            9912399*ones(size(outMat1(:, 1))),...							% blank
            (outMat1(:, 8)./ repmat(outMat1(:, 1), 1, 1) *100000) ,...		% TBinc
            9912399*ones(size(outMat1(:, 1))),...							% blank
            (outMat1(:, 14)./ repmat(outMat1(:, 1), 1, 1) *100000) ,...     % TBmortExcess
            9912399*ones(size(outMat1(:, 1))),...							% blank
            outMat1(:, 11),...                                              % numPplCAT12
            9912399*ones(size(outMat1(:, 1))),...							% blank
            outMat1(:, 12),...                                              % numPplCat4
            outMat1(:, 9:10),...                                            % fracIncIsMDR, fracTBincUnder15
            outMat2(:, 6) ./ outMat1(:, 1),...                              % proportion of pop TBlatPrev
            9912399*ones(size(outMat1(:, 1))),...							% not doing prop disease due to recent infecton
            outMat1(:, 13)./ outMat1(:, 1),...                              % annRiskofInfection
            9912399*ones(size(outMat1(:, 1))),...							% not doing R0
            (outMat2(:, 2:3)./ repmat(outMat1(:, 1), 1, 2) *100000) ...     % TBmdrPrev,TBprev
            ];
        header = 'Year,  pop, blank, totalNotifications,  blank, TBinc, blank, TBmort_everTreatment, blank, numPplCAT12, blank, numPplCat4,  fracIncIsMDR, fracTBincUnder15,fracIsLatPrev, blank, annRiskofInfection, TBmdrPrev,TBprev';
        cd(masterFolder)
        tablePrinter(header, formattedOutMatFinal, 'TB_MAC_Output_formatted', currentFold);
        
        %make a vector of the appropriate population
        actualPopInTime = TBparams.IndianPop(find(ismember(popYear,[1996:1:2035])))';
        
        %make epi results
        TB_MAC_epi_results = 12345678*ones(40,11);
        TB_MAC_epi_results(:,1) =  outMat1(:, [18])./ outMat1(:, 1) *100000;% TB_notif   TB notification (n/100k/year) (adults, all types, all sectors)
        TB_MAC_epi_results(:,2) =  outMat2(:, 9)   ./ outMat1(:, 1) *100000; %         TB_prev  TB disease prevalence (n/100k/) (adults, all types)
        TB_MAC_epi_results(:,3) =  ((1-outMat1(:, 15)) .* (outMat1(:, 6))./ outMat1(:, 1) *100000); %         TB_mort.   TB mortality (n/100k/year) (adults, all types, all sectors)
        TB_MAC_epi_results(:,4) =  (1-outMat1(:, 10)) .* (outMat1(:, 8)./ outMat1(:, 1) *100000) ; %         TB_incid  TB disease incidence (n/100k/year)  (adults, all types)
        TB_MAC_epi_results(:,5) =  (1-outMat2(:,7)) .* actualPopInTime; %         Pop_size_adults
        TB_MAC_epi_results(:,6) =  outMat2(:,10); %         MDR_per_new_incid  % of MDR TB in incident new and retreatment cases
        TB_MAC_epi_results(:,7) =  outMat2(:,11); %         MDR_per_retreat_incid
        TB_MAC_epi_results(:,8) =  outMat2(:, 6) ./ outMat1(:, 1); %         Prop_LTBI in total pop
        TB_MAC_epi_results(:,9) =  outMat2(:,12); %         Prop_TB_recent   Proportion of TB disease due to recent infection (w/in 2 yrs)
        TB_MAC_epi_results(:,10) = outMat1(:, 13)./ outMat1(:, 1) ; %         ARI  (n new infections/total population)
        TB_MAC_epi_results(:,11) = outMat2(:,14); %         Avg_time_diag
        header = 'TB_notif,TB_prev,TB_mort,TB_incid,Pop_size_adults,MDR_per_new_incid,MDR_per_retreat_incid,Prop_LTBI,Prop_TB_recent,ARI,Avg_time_diag';
        tablePrinter(header, TB_MAC_epi_results, 'TB_MAC_epi_results', currentFold);
        
        %make Daly1, 2, 3 results
        row_yr2034 = 39;
        TB_MAC_daly1_results = 12345678*ones(40,2); TB_MAC_daly2_results  = 12345678*ones(40,10);  TB_MAC_daly3_results= 12345678*ones(2,10);
        TB_MAC_daly1_results(:,1) = (outMat2(:, 3)./ outMat1(:, 1)) .* actualPopInTime; %Population totals at mid-year, by health state Active TB	No Active TB
        TB_MAC_daly1_results(:,2) = actualPopInTime - TB_MAC_daly1_results(:,1); %Population totals at mid-year, by health state Active TB	No Active TB
        TB_MAC_daly2_results = outMat4(:,1:10) ./ repmat(outMat1(:, 1), 1, 10) .* repmat(actualPopInTime, 1, 10);% Table 2: Total deaths in each year, by age group Age 0-9	Age 10-19	Age 20-29	Age 30-39	Age 40-49	Age 50-59	Age 60-69	Age 70-79	Age 80-89	Age 90+
        TB_MAC_daly3_results = [outMat4(row_yr2034, 11:20) ./ repmat(outMat1(row_yr2034, 1), 1, 10) .* repmat(TBparams.IndianPop(popYear == 2034), 1, 10);...
            outMat4(row_yr2034, 21:30)/ repmat(outMat1(row_yr2034, 1), 1, 10) .* repmat(TBparams.IndianPop(popYear == 2034), 1, 10)];% Table 3: Population at end 2034 disaggregated by age and health state
        tablePrinter('Pop totals Active TB,	Pop No Active TB', TB_MAC_daly1_results, 'TB_MAC_daly1_results', currentFold);
        tablePrinter('Total deaths Age 0-9,Age 10-19,Age 20-29,Age 30-39,Age 40-49,Age 50-59,Age 60-69,Age 70-79,Age 80-89,Age 90+', TB_MAC_daly2_results, 'TB_MAC_daly2_results', currentFold);
        tablePrinter('total alive in 2034 by age, first row is actTBalive_ages, second row is noTB', TB_MAC_daly3_results, 'TB_MAC_daly3_results', currentFold);
        
        %make econ results
        TB_MAC_econ_results = 12345678*ones(40,20);
        TB_MAC_econ_results(:,3) = outMat1(:,19)        ./ outMat1(:, 1) .* actualPopInTime; %1a+b. Total person months with untreated Active TB
        TB_MAC_econ_results(:,5) = sum(outMat1(:,2:3),2)./ outMat1(:, 1) .* actualPopInTime ; %Total TB diagnoses
        TB_MAC_econ_results(:,6) =  outMat1(:,20)       ./ outMat1(:, 1) .* actualPopInTime;% 3a. Total suspects screened via ACF
        TB_MAC_econ_results(:,7) = outMat1(:,21)        ./ outMat1(:, 1) .* actualPopInTime ;
        TB_MAC_econ_results(:,8) = outMat1(:,4)         ./ outMat1(:, 1) .* actualPopInTime ;
        TB_MAC_econ_results(:,9) = (outMat1(:,24) + outMat1(:,29))  ./ outMat1(:, 1) .* actualPopInTime ; %person months in CatI + private trtNaive
        TB_MAC_econ_results(:,10) = (outMat1(:,25) + outMat1(:,28)) ./ outMat1(:, 1) .* actualPopInTime ; %person months in CatII+ private trtExperienced
        TB_MAC_econ_results(:,11) = TB_MAC_econ_results(:,9) + TB_MAC_econ_results(:,10);
        TB_MAC_econ_results(:,12) = outMat1(:,26)       ./ outMat1(:, 1) .* actualPopInTime ; %person months in CatIV
        TB_MAC_econ_results(:,13) = outMat1(:,27)       ./ outMat1(:, 1) .* actualPopInTime ; %person months in private
        TB_MAC_econ_results(:,14) = 0;  %no MDR in private
        
        TB_MAC_econ_results(:,15) = outMat1(:,22)       ./ outMat1(:, 1) .* actualPopInTime ;
        TB_MAC_econ_results(:,16) = outMat1(:,23)       ./ outMat1(:, 1) .* actualPopInTime ;
        header = 'TB_mo_pre_diag_atmpt,TB_mo_post_diag_atmpt,TB_mo_total,TB_suspects_passive,TB_diagnoses_passive,TB_suspects_ACF,TB_diagnoses_ACF,DST,1st_line_tx_naiv,1st_line_tx_expd,1st_line_all,MDR_reg,fract_1st_line_priv,fract_MDR_reg_priv,LTBI_screen,LTBI_treatment,IPT_screen,IPT_treatment,ART_init,ART';
        tablePrinter(header, TB_MAC_econ_results, 'TB_MAC_econ_results', currentFold);
        
        %make comparison graphs
        if size(TB_MACtargetFolds,2) > 1
            toPlot{1} = [toPlot{1}, formattedOutMatFinal(:, 4)];  %notifications
            toPlot{2} = [toPlot{2}, formattedOutMatFinal(:, 6)];  %TBinc
            toPlot{3} = [toPlot{3}, formattedOutMatFinal(:, 8)];  %TBmortExcess
            toPlot{4} = [toPlot{4}, formattedOutMatFinal(:, 10)];  %numPplCAT12
            toPlot{5} = [toPlot{5}, formattedOutMatFinal(:, 12)];  %numPplCat4
        end
        
        %%%FOR THE OUTPUT TABLES
        if TB_MACtargets_compiled == 1    %put all the interventions into one big table
            
            %get the correct intervention Num and scenarioStrengthStr
            % For column A (Intervention), please use the following codes: 0 (i.e. the basecase), 1, 2, 2a, 2b, 2c, 3a, 3b, 4, 5, 6, and 7.
            % For column B (Scenario), please use the following codes: base_case, cty_25, cty_50, cty_75, cty_100 (i.e. different levels of the country expert intervention description), and adv (the advocate intervention description)
            yearRange = [21:40]'; %or 20:40 for base case.  Year 1 is 1996.
            scenarioFolderName  = TB_MACtargetFolds{folderNum};
            scenarioFolderName = scenarioFolderName(1:end-1);  %take off the slash at the end
            scaleNum = scenarioFolderName(end-3: end);  %p25, 0p5, p75, or letters (base,etc)
            if (size(scenarioFolderName, 2) == 4)  %is base case
                interventionNum = '0';  %
                ScenarioStrengthStr = 'base_case';
                yearRange = [5:40]';
            else
                interventionNum = scenarioFolderName(1: 7);  %TBMAC1_, etc or TBMAC3b, etc
                interventionNum = strrep(strrep(interventionNum, '_', ''), 'TBMac', '');  %leave only the intervention number
                if strcmpi(interventionNum, '6') == 1  %my scenario 6_combination is actually intervention 7
                    interventionNum = '7';
                elseif strcmpi(interventionNum, '3') == 1
                    interventionNum = '3a';
                end
                ScenarioStrengthStr = 'cty_100';
                if strcmpi(scaleNum, '0p25') == 1
                    ScenarioStrengthStr = 'cty_25';
                elseif strcmpi(scaleNum, '_0p5') == 1
                    ScenarioStrengthStr = 'cty_50';
                elseif strcmpi(scaleNum, '0p75') == 1
                    ScenarioStrengthStr = 'cty_75';
                elseif strcmpi(scaleNum, '_adv') == 1
                    ScenarioStrengthStr = 'adv';
                end
            end

            CountryStr = 'india';
            dataLength = size(yearRange,1);  % For column D (Year), please use years 2015 through 2035 inclusive for the basecase, and years 2016 through 2035 for all other interventions and scenarios.
            
            colA = repmat({interventionNum},dataLength,1);
            colB = repmat({ScenarioStrengthStr},dataLength,1);
            colC = repmat({CountryStr},dataLength,1);
            yearVals = [1996:1:2036]';
            colD = yearVals(yearRange);
            colD_ages = {'0_9';'10_19';'20_29';'30_39';'40_49';'50_59';'60_69';'70_79';'80_89';'90+'};
            
            % TB_MAC_epi_results
            T1 = table(colA,colB,colC,colD,TB_MAC_epi_results(yearRange,:)); TB_MAC_epi_results_all = [TB_MAC_epi_results_all;T1];

            % TB_MAC_daly1_results
            T1 = table(colA,colB,colC,colD,TB_MAC_daly1_results(yearRange,:)); TB_MAC_daly1_results_all = [TB_MAC_daly1_results_all;T1];
            
            % TB_MAC_daly2_results
            T1 = table(colA,colB,colC,colD,TB_MAC_daly2_results(yearRange,:)); TB_MAC_daly2_results_all = [TB_MAC_daly2_results_all;T1];
            
            % TB_MAC_daly3_results
            T1 = table(colA(1:10),colB(1:10),colC(1:10),colD_ages,TB_MAC_daly3_results'); TB_MAC_daly3_results_all = [TB_MAC_daly3_results_all;T1];
            
            % TB_MAC_econ_results
            T1 = table(colA,colB,colC,colD,TB_MAC_econ_results(yearRange,:)); TB_MAC_econ_results_all = [TB_MAC_econ_results_all;T1];
        end
    end  %end folderNum loop
    
    if TB_MACtargets_compiled == 1  %write the results tables out
        cd(fullfile(masterFolder,'paperFigures'))
        writetable(TB_MAC_epi_results_all,'TB_MAC_epi_results_all.csv');
        writetable(TB_MAC_daly1_results_all,'TB_MAC_daly1_results_all.csv');
        writetable(TB_MAC_daly2_results_all,'TB_MAC_daly2_results_all.csv');
        writetable(TB_MAC_daly3_results_all,'TB_MAC_daly3_results_all.csv');
        writetable(TB_MAC_econ_results_all,'TB_MAC_econ_results_all.csv');
    end
    
    if size(TB_MACtargetFolds,2) > 1   %make a graph for easy comparison
        titleArray = {'Notifications (n out of 100,000)', 'Incidence (n out of 100,000)', 'Mortality (n out of 100,000)','Person Months on DOTS','Person Months on MDR-TB Treatment'};
        for TBMACoutputMeas = 1:5
            subplot(3,2,TBMACoutputMeas)
            plot([2000:1:2035]', toPlot{TBMACoutputMeas}(5:end,:)); %1:4))
            xlabel('Year');
            title(titleArray{TBMACoutputMeas});
            if TBMACoutputMeas == 1
                legend('base','Increase DOTS', 'Improve DOTS', 'Xpert', 'ACF', 'ACF w/ Latent Tx','Combination', 'Location', 'NorthWest');
            end
        end
        titleStr = 'allTBmacInterventions';
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 25])
        print(printFormat,'-r300',[outputFolder '/' titleStr]);
    end
    
end



