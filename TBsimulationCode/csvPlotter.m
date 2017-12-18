function csvPlotter(csvFile, folderName, plotResolution)

%function csvPlotter(csvFile, folderName, plotResolution)

%

%csvFile is the output of csvPlotMatrix, is a string

%folderName is where you want the csvFiles to print

%plotResolution should be '-r70'





currentDirectory = pwd;

cd(folderName);  %go to the output directory

fullPlottingMatrix = dlmread(csvFile, ',' ,0,0);

% disp(sprintf('plotting matrix is %i by %i',size(fullPlottingMatrix,1),size(fullPlottingMatrix,2)    ));

% size(fullPlottingMatrix)

cd(currentDirectory); %get back out to the original directory



%make indicies easier to read; also note that dead and not yet born will be wrong if recno is not decompressed

numHealthOutcomes = 14;

healthy = 1; latentSens = 2; latentMdr = 3; activeSens =4; activeMdr = 5; notYetBorn = 6; dead = 7;

catI_sens = 8; catII_sens = 9; catIII = 10; catIV_sens = 11; catI_MDR = 12; catII_MDR = 13; catIV_MDR = 14;



male = 1; female = 2; 

smoking = 1; nonSmoking = 2; 

ageBinMat_1 = 1; ageBinMat_2 = 2; ageBinMat_3 = 3;



numCol = size( fullPlottingMatrix,2);



%--------------%%GRAB CHARACTERISTIC COLUMNS%%-------------------%



sexTypes = 2;

maleSelector = [   true(1, numCol/sexTypes )   ,   false(1, numCol/sexTypes )   ];

femaleSelector = (~maleSelector);



smokingTypes = 2;

smokingDivisor = sexTypes*smokingTypes;

oneSetSmokers =[ true(1, numCol/smokingDivisor ) , false( 1, numCol/smokingDivisor) ];

smokingSelector = [ oneSetSmokers, oneSetSmokers  ];

nonSmokingSelector = (~smokingSelector);



ageTypes = 3;

ageDivisor = smokingDivisor * ageTypes;

oneSetAge = [true(1, numCol/ageDivisor ) , false(1, numCol/ageDivisor ),  false(1, numCol/ageDivisor )];

age1_selector = logical([oneSetAge, oneSetAge, oneSetAge, oneSetAge]);

oneSetAge = [false(1, numCol/ageDivisor ) , true(1, numCol/ageDivisor ),  false(1, numCol/ageDivisor )];

age2_selector = [oneSetAge, oneSetAge, oneSetAge, oneSetAge];

age3_selector = (~age1_selector & ~age2_selector);



%healthy is every 14

nothingSelector = false(1,numCol);

selectorName = {'healthy','latentSens'};

for  i = 1:14

    selectorName{i} = nothingSelector;

    selectorName{i}(i:14: numCol ) = 1;

end



healthySelector = selectorName{healthy};

latentSensSelector = selectorName{latentSens};

latentMdrSelector = selectorName{latentMdr};

activeSensSelector = selectorName{activeSens};

activeMdrSelector = selectorName{activeMdr};

notYetBornSelector = selectorName{notYetBorn};

deadSelector = selectorName{dead};

catI_sensSelector = selectorName{catI_sens};

catII_sensSelector = selectorName{catII_sens};

catIIISelector = selectorName{catIII};

catIV_sensSelector = selectorName{catIV_sens};

catI_MDRSelector = selectorName{catI_MDR};

catII_MDRSelector = selectorName{catII_MDR};

catIV_MDRSelector = selectorName{catIV_MDR};



aliveSelector = (healthySelector | latentSensSelector | latentMdrSelector |  ...
    activeSensSelector | activeMdrSelector | catI_sensSelector | ...
    catII_sensSelector | catIIISelector | catIV_sensSelector | catI_MDRSelector | ...
    catII_MDRSelector | catIV_MDRSelector);



%--------------%%MAKE GRAPHS%%-------------------%



for subsetGrp = 1:20



    if subsetGrp == 1

        subset = aliveSelector;

        titleSuffix = '';        

    elseif subsetGrp == 2

        subset = (maleSelector & aliveSelector);

        titleSuffix = '_maleOnly';

    elseif subsetGrp == 3

        subset = (femaleSelector & aliveSelector);

        titleSuffix = '_femaleOnly';

    elseif subsetGrp == 4

        subset = (smokingSelector & aliveSelector);

        titleSuffix = '_smokingOnly';

    elseif subsetGrp == 5

        subset = (nonSmokingSelector & aliveSelector);

        titleSuffix = '_nonSmokingOnly';

    elseif subsetGrp == 6

        subset = (age1_selector & aliveSelector);

        titleSuffix = '_age0_20Only';

    elseif subsetGrp == 7

        subset = (age2_selector & aliveSelector);

        titleSuffix = '_age21_60Only';

    elseif subsetGrp == 8

        subset = (age3_selector & aliveSelector);

        titleSuffix = '_age61_100Only';

        

    elseif subsetGrp == 9

        subset = (age1_selector & smokingSelector & aliveSelector);

        titleSuffix = '_smoking_age0_20Only';

    elseif subsetGrp == 10

        subset = (age2_selector & smokingSelector & aliveSelector);

        titleSuffix = '_smoking_age21_60Only';

    elseif subsetGrp == 11

        subset = (age3_selector & smokingSelector & aliveSelector);

        titleSuffix = '_smoking_age61_100Only';

    elseif subsetGrp == 12

        subset = (age1_selector & nonSmokingSelector & aliveSelector);

        titleSuffix = '_nonSmoking_age0_20Only';

    elseif subsetGrp == 13

        subset = (age2_selector & nonSmokingSelector & aliveSelector);

        titleSuffix = '_nonSmoking_age21_60Only';

    elseif subsetGrp == 14

        subset = (age3_selector & nonSmokingSelector & aliveSelector);

        titleSuffix = '_nonSmoking_age61_100Only';

    

    elseif subsetGrp == 15

        subset = (age1_selector & maleSelector & aliveSelector);

        titleSuffix = '_male_age0_20Only';

    elseif subsetGrp == 16

        subset = (age2_selector & maleSelector & aliveSelector);

        titleSuffix = '_male_age21_60Only';

    elseif subsetGrp == 17

        subset = (age3_selector & maleSelector & aliveSelector);

        titleSuffix = '_male_age61_100Only';

    elseif subsetGrp == 18

        subset = (age1_selector & femaleSelector & aliveSelector);

        titleSuffix = '_female_age0_20Only';

    elseif subsetGrp == 19

        subset = (age2_selector & femaleSelector & aliveSelector);

        titleSuffix = '_female_age21_60Only';

    elseif subsetGrp == 20

        subset = (age3_selector & femaleSelector & aliveSelector);

        titleSuffix = '_female_age61_100Only';

        

    end



    %%%%%%%%%%%%%health outcomes graph%%%%%%%%%%

    aliveTots = sum(fullPlottingMatrix(:,(subset & aliveSelector)),2);

    healthyTots = sum(fullPlottingMatrix(:,(subset & healthySelector)),2);

    latSensTots = sum(fullPlottingMatrix(:,(subset & latentSensSelector)),2);

    latMDRTots = sum(fullPlottingMatrix(:,(subset & latentMdrSelector)),2);

    activeSensTots = sum(fullPlottingMatrix(:,(subset & (activeSensSelector | catI_sensSelector | catII_sensSelector | catIV_sensSelector ))),2);

    activeMdrTots = sum(fullPlottingMatrix(:,(subset & (activeMdrSelector | catI_MDRSelector | catII_MDRSelector |  catIV_MDRSelector))),2);



    toPlot = [healthyTots, latSensTots, latMDRTots, activeSensTots, activeMdrTots];

    legendText = {'healthy','latent DS','latent MDR','active DS', 'active MDR'};

    propMat = toPlot ./ [aliveTots, aliveTots, aliveTots, aliveTots, aliveTots];

    titleSt = strcat('HealthOutcomes', titleSuffix);

    tableTitle = strcat(titleSt, '.png');

    tableHeader = 'healthy,latentSens,latentMdr,actSens, actMdr, prop_healthy, prop_latentSens, prop_latentMdr, prop_actSens, prop_actMDR';

    tableMat = [toPlot,propMat];



    %PLOT

    subplot(2,1,1);

    bar(toPlot,'stack')

    leg = legend(legendText);

    set(leg,'Location','NorthEastOutside');

    title(titleSt);

    xlabel('Year')

    ylabel('No. of People')

    hold on;

    subplot(2,1,2);

    bar(propMat,'stack');

    ylim([0,1]);

    title('(by proportion of subset, if there is one)')

    xlabel('Year');

    ylabel('Proportion of subset');

    print('-dpng',plotResolution,[folderName '/' tableTitle]);

    tablePrinter(tableHeader, tableMat, titleSt, folderName);

    close all



    %only plot post treatment

    
    x = [1966:5:2046];
    toPlot = [healthyTots(20:end,:), latSensTots(20:end,:), latMDRTots(20:end,:), activeSensTots(20:end,:), activeMdrTots(20:end,:)];
    
    x = x(1:size(toPlot,1));
    legendText = {'healthy','latent DS','latent MDR','active DS', 'active MDR'};

    

    propMat = toPlot ./ [aliveTots(20:end,:), aliveTots(20:end,:), aliveTots(20:end,:), aliveTots(20:end,:), aliveTots(20:end,:)];

    titleSt = strcat('HealthOutcomes_postBurnIn', titleSuffix);

    tableTitle = strcat(titleSt, '.png');

    tableHeader = 'healthy,latentSens,latentMdr,actSens, actMdr, prop_healthy, prop_latentSens, prop_latentMdr, prop_actSens, prop_actMDR';

    tableMat = [toPlot,propMat];

    %PLOT

    subplot(2,1,1);

    bar(x,toPlot,'stack')

    leg = legend(legendText);

    set(leg,'Location','NorthEastOutside');

    title(titleSt);

    xlabel('Year')

    ylabel('No. of People')

    hold on;

    subplot(2,1,2);

    bar(x,propMat,'stack');

    ylim([0,1]);

    title('(by proportion of subset, if there is one)')

    xlabel('Year');

    ylabel('Proportion of Subset');

    print('-dpng',plotResolution,[folderName '/' tableTitle]);

    tablePrinter(tableHeader, tableMat, titleSt, folderName);

    close all

    

    

    if subsetGrp == 1

        %also plot prevalence graphs with everyone (no subsetting)

        toPlot = [activeSensTots(20:end,:), activeMdrTots(20:end,:)] ./ [aliveTots(20:end,:), aliveTots(20:end,:)];

        bar(x,toPlot, 'stack');

        legend('actSensPrevalence', 'actMDRprevalence');

        title('TB prevalence');

        xlabel('Year');

        ylabel('Prevalence');

        print('-dpng',plotResolution,[folderName '/' 'TBprevalence.png']);

        tablePrinter('actSensPrevalence, actMDRPrevalence', toPlot, 'TBprevalence', folderName);

        close all

        

        %do prevalence by rate (ie, total active out of 100,000) and put

        %the WHO numbers on too for comparison

        dataYear = [1990	; 1991	; 1992	; 1993	; 1994	;  1995	; 1996	; 1997	; 1998	; 1999	; 2000	; 2001	; 2002	; 2003	; 2004	; 2005	; 2006	; 2007	; 2008	; 2009	; 2010	];

        actTBPrev = [459 ; 460 ; 460	; 461	; 462	; 462	; 463	; 464	; 464	; 465	; 466	; 466	; 436	; 409	; 383	; 358	; 335	; 314	; 294	; 275	; 256];

        distToUpperLimit = [56	; 56	; 56	; 56	; 56	; 57	; 56	; 56	; 57	; 56	; 56	; 57	; 62	; 66	; 71	; 78	; 84	; 91	; 99	; 107	; 117];

        distToLowerLimit = [52	; 53	; 52	; 53	; 53	; 53	; 53	; 53	; 53	; 53	; 54	; 53	; 57	; 62	; 66	; 70	; 74	; 80	; 85	; 90	; 95];   

        

        %make simulation prevalences into rates out of 100000

        smallx = [1991:5:2011];

        toPlot = ((activeSensTots(20:end,:)+activeMdrTots(20:end,:)) ./ aliveTots(20:end,:)) *100000;

        %p = plot(x,toPlot);

        p = plot(smallx,toPlot(6:10));

        set(p,'Color','red','LineWidth',1);

        title('Total active TB, out of 100000, with WHO estimates');

        xlabel('Year');

        ylabel('Prevalence Rate');

        hold on;

        errorbar(dataYear, actTBPrev,distToLowerLimit,distToUpperLimit);

        ylim([0 600]);

        print('-dpng',plotResolution,[folderName '/' 'TBprevWHOcomparison.png']);

        close all

        

        %do mdr as a percentage of total active and put WHO numbers on

        toPlot = activeMdrTots(20:end,:) ./ (activeMdrTots(20:end,:) + activeSensTots(20:end,:));

        %p = plot(x,toPlot);

        p = plot(smallx,toPlot(6:10));

        set(p,'Color','red','LineWidth',1);

        title('Percentage MDR TB Out of Total Active TB, with WHO estimates');

        hold on;

        %data from WHO MDR Report 2004 - 2007

        mdrPrevYear = [1997; 1999; 1999; 2001;2001; 2001;2004; 2006];  %these are calculated with the numbers above on WHO_mdr_fracOfTotalPop in WHO articles fold (numMDR/numSusceptible)

        mdrPercOfActTB = [4.166666667;3.225806452;3.921568627;0.74906367;0.632911392;3.652968037;2.727272727;2.993527508];  

        scatter(mdrPrevYear,0.01*mdrPercOfActTB, 'o');

        xlabel('Year');

        ylabel('Proportion MDR Out of Total Active TB');        

        print('-dpng',plotResolution,[folderName '/' 'MDRpercWHOcomparison.png']);

        close all

        

        %do cases of mdr as a percentage of total Indian population put WHO numbers on

        toPlot = activeMdrTots(20:end,:) ./ sum([healthyTots(20:end,:), latSensTots(20:end,:), latMDRTots(20:end,:), activeSensTots(20:end,:), activeMdrTots(20:end,:)], 2);

        p = plot(x,toPlot);

        set(p,'Color','red','LineWidth',1);

        title('Percentage MDR TB Out of Total Population, with WHO estimates');

        hold on;

        %data from WHO MDR Report 2004 - 2007

        mdrPrevYear = [1997; 1999; 1999; 2001;2001; 2001;2004; 2006];  %these are calculated with the numbers above on WHO_mdr_fracOfTotalPop in WHO articles fold

        mdrPercOfTotal = [0.000193333;0.00015;0.000182353;0.000034906;0.000029494;0.000170228;0.000104455;0.000100283];



        scatter(mdrPrevYear,mdrPercOfTotal, 'o');

        xlabel('Year');

        ylabel('Proportion MDR Out of Total Active TB');

        print('-dpng',plotResolution,[folderName '/' 'MDRpercWHOcomparison_ofTotPop.png']);

        tablePrinter('activeMdrTots', toPlot, 'MDR_TBprevalence', folderName);

        close all

        

    end



    %%%%%%%%%%%%%activeTB graph%%%%%%%%%%

    activeSensTots = sum(fullPlottingMatrix(:,(subset & (activeSensSelector))),2);

    activeSensCatITots = sum(fullPlottingMatrix(:,(subset & (catI_sensSelector))),2);

    activeSensCatIITots = sum(fullPlottingMatrix(:,(subset & (catII_sensSelector))),2);

    activeSensCatIVTots = sum(fullPlottingMatrix(:,(subset & (catIV_sensSelector))),2);

    activeMdrTots = sum(fullPlottingMatrix(:,(subset & (activeMdrSelector))),2);

    activeMdrCatITots = sum(fullPlottingMatrix(:,(subset & ( catI_MDRSelector))),2);

    activeMdrCatIITots = sum(fullPlottingMatrix(:,(subset & ( catII_MDRSelector))),2);

    activeMdrCatIVTots = sum(fullPlottingMatrix(:,(subset & ( catIV_MDRSelector))),2);



    toPlot = [activeSensTots(20:end,:), activeSensCatITots(20:end,:), activeSensCatIITots(20:end,:), activeSensCatIVTots(20:end,:), activeMdrTots(20:end,:), activeMdrCatITots(20:end,:), activeMdrCatIITots(20:end,:), activeMdrCatIVTots(20:end,:)];

    actTBTots = sum(toPlot,2);

    legendText = {'DS untreated', 'DS CatI', 'DS CatII', 'DS CatIV', 'MDR untreated', 'MDR CatI', 'MDR CatII', 'MDR CatIV'};

    propMat = toPlot ./ [actTBTots, actTBTots, actTBTots, actTBTots, actTBTots, actTBTots, actTBTots, actTBTots];    

    titleSt = strcat('actTBoutcomes', titleSuffix);

    tableTitle = strcat(titleSt, '.png');

    tableHeader = 'actSens, actSensCatI, actSensCatII, actSensCatIV, actMdr, actMdrCatI, actMdrCatII, actMdrCatIV, prop_actSens, prop_actSensCatI, prop_actSensCatII, prop_actSensCatIV, prop_actMdr, prop_actMdrCatI, prop_actMdrCatII, prop_actMdrCatIV';

    tableMat = [toPlot,propMat];



    %PLOT

    subplot(2,1,1);

    bar(x,toPlot,'stack')

    leg = legend(legendText);

    set(leg,'Location','NorthEastOutside');

    title(titleSt);

    xlabel('Years');

    ylabel('Number of People with Active TB');

    hold on;

    subplot(2,1,2);

    bar(x,propMat,'stack');

    ylim([0,1]);

    title('(by proportion of subset with active TB)')

    xlabel('Years');

    ylabel('Proportion of subset with Active TB');

    print('-dpng',plotResolution,[folderName '/' tableTitle]);

    tablePrinter(tableHeader, tableMat, titleSt, folderName);

    close all



end











