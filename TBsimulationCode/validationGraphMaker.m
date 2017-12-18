%model validation graphs
masterFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\May28\';

%WHAT DO YOU WANT TO RUN????
populationGraph = 1;  %turn population growth graph on and off
if populationGraph == 1
	%needed for the populationGrowthGraph
	baseCase = 'base\p01_2014-05-28_16-10-36';
end

lifeExpectancyCIs = 0;  %turn the life expectancy tables on and off
if lifeExpectancyCIs == 1
	%needed for the lifeExpectancyCIs
	calibration1990 = 'noTreat_fullCatIV_empUpta_calibration1990\';
	fold1990 = {'Cal_2013-02-13_11-29-11';'Cal_2013-02-13_11-33-48';'Cal_2013-02-13_11-38-27'};
	n1990 = size(fold1990,1);
	calibration2000 = 'noTreat_fullCatIV_empUpta_calibration2000\';
	fold2000 = {'Cal_2013-02-13_11-43-08';'Cal_2013-02-13_11-47-50';'Cal_2013-02-13_11-52-31'};
	n2000 = size(fold2000,1);
end

treatmentStatistics = 1;  %turn on the treatment statistics maker
if treatmentStatistics == 1
	%needed for the treatmentStatistics
	baseCase = 'base\p01_2014-05-28_16-10-36';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%make comparisons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotResolution = '-r70';
plotFixerFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs'; %this is where plotfixer lives;
%validationGraphsFold = strcat(masterFold, 'validationGraphs');  %where the graphs should be outputted to
validationGraphsFold = strcat(masterFold, 'paperFigures');  %where the graphs should be outputted to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if populationGraph == 1
    
    %grab data
    cd(strcat(masterFold, baseCase))
    basecase_ageStructure = dlmread('AgeStructure.csv', ',' ,1,0);
    
    
    %POPULATION GROWTH
    WHO_popGrowthRate = [...
        1995	0.020760485
        2000	0.01854086
        2005	0.016347882
        2010	0.014836458
        2015	0.013654425
        2020	0.012029772
        2025	0.010389867
        2030	0.008845217
        2035	0.007393589
        2040	0.005978851
        ];
    
    WHO_pop_inThou = [...
        1995	964486
        2000	1053898
        2005	1140043
        2010	1224614
        2015	1308221
        2020	1386909
        2025	1458958
        2030	1523482
        2035	1579802
        2040	1627029
        ];
    
    simulationPop = sum(basecase_ageStructure(26:end-1,1:5),2);
    
    for i = 2: size(simulationPop,1)
        WHO_fiveYrGrowthRate(i) = (WHO_pop_inThou(i,2)-WHO_pop_inThou(i-1,2))/(WHO_pop_inThou(i-1,2));
        simu_fiveYrGrowthRate(i) = (simulationPop(i,1)-simulationPop(i-1,1))/(simulationPop(i-1,1));
    end
    
    plot(WHO_pop_inThou(1:size(simulationPop,1),1), WHO_fiveYrGrowthRate(1:end-1), 'bo');
    hold on;
    plot(WHO_pop_inThou(1:size(simulationPop,1),1), simu_fiveYrGrowthRate, 'r.');
    ylim([0 0.5]);
    legend('UN Population Projection Data', 'Simulation Data')
    cd(plotFixerFold);
    plotfixer;
    print('-dpng',plotResolution,[validationGraphsFold '/'  'PopulationGrowthValidation']);
    close all
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lifeExpectancyCIs == 1
    
    %make 1990 and 2000 mean and 95% CI for
    
    %%%run ShowFishes
    TBactivated_1990 = zeros(n1990, 5);
    femaleLE_1990 = zeros(n1990, 5);
    maleLE_1990 = zeros(n1990, 5);
    censored_1990 = zeros(n1990, 2);
    TBactivated_2000 = zeros(n2000, 5);
    femaleLE_2000= zeros(n2000, 5);
    maleLE_2000 = zeros(n2000, 5);
    censored_2000= zeros(n2000, 2);
    
    for i = 1:n1990
        cd(masterFold)
        folder = strcat(calibration1990,fold1990{i});
        ShowFishes_mortalityDebugger(folder, 1990);
        cd(folder)
        postTB = dlmread('postActivationLifeExp.csv', ',' ,1,0);  %postTB LE, nonsmokers (given live to age30), smokers(given live to age30), percActivate, totNumPpl
        lifeExp = dlmread('LifeExp.csv', ',' ,1,0);  %WHO LE	 Simu LE	 WHO LE from 30	 nonsmokers Simu LE from 30	 smokers Simu LE from 30	 row1 male and row2 female
        nonZeroMonSinceActTruncations = dlmread('nonZeroMonSinceActTruncations.csv', ',' ,1,0); %number of nonzero elements in MonSinceActTruncations
        
        TBactivated_1990(i,:) = postTB;
        femaleLE_1990(i,:) = lifeExp(2,:);
        maleLE_1990(i,:) = lifeExp(1,:);
        censored_1990(i,:) = nonZeroMonSinceActTruncations';
    end
    
    for i = 1:n2000
        cd(masterFold)
        folder = strcat(calibration2000,fold2000{i});
        ShowFishes_mortalityDebugger(folder, 2000);
        cd(folder)
        postTB = dlmread('postActivationLifeExp.csv', ',' ,1,0);  %postTB LE, nonsmokers (given live to age30), smokers(given live to age30), percActivate, totNumPpl
        lifeExp = dlmread('LifeExp.csv', ',' ,1,0);  %WHO LE	 Simu LE	 WHO LE from 30	 nonsmokers Simu LE from 30	 smokers Simu LE from 30	 row1 male and row2 female
        nonZeroMonSinceActTruncations = dlmread('nonZeroMonSinceActTruncations.csv', ',' ,1,0); %number of nonzero elements in MonSinceActTruncations
        
        TBactivated_2000(i,:) = postTB;
        femaleLE_2000(i,:) = lifeExp(2,:);
        maleLE_2000(i,:) = lifeExp(1,:);
        censored_2000(i,:) = nonZeroMonSinceActTruncations;
    end
    
    
    disp('these are the 1990 censored')
    censored_1990
    
    disp('these are the 2000 censored')
    censored_2000
    
    cd(masterFold)
    
    TBactivated_1990_means =  mean(TBactivated_1990,1);
    [a,b,TBactivated_1990_CIs] = ttest(TBactivated_1990);
    femaleLE_1990_means = mean(femaleLE_1990,1);
    [a,b,femaleLE_1990_CIs] = ttest(femaleLE_1990);
    maleLE_1990_means = mean(maleLE_1990,1);
    [a,b,maleLE_1990_CIs] = ttest(maleLE_1990);
    
    TBactivated_2000_means =  mean(TBactivated_2000,1);
    [a,b,TBactivated_2000_CIs] = ttest(TBactivated_2000);
    femaleLE_2000_means =  mean(femaleLE_2000,1);
    [a,b,femaleLE_2000_CIs] = ttest(femaleLE_2000);
    maleLE_2000_means = mean(maleLE_2000,1);
    [a,b,maleLE_2000_CIs] = ttest(maleLE_2000);
    
    meansAndCIs_1990 = [ TBactivated_1990_means , maleLE_1990_means , femaleLE_1990_means ; TBactivated_1990_CIs , maleLE_1990_CIs, femaleLE_1990_CIs ];
    meansAndCIs_2000 = [ TBactivated_2000_means , maleLE_2000_means , femaleLE_2000_means ; TBactivated_2000_CIs , maleLE_2000_CIs, femaleLE_2000_CIs ];
    
    rowLabelStr = 'postTB_LE, postTB_LEnonsmokers_postAge30, postTB_LEsmokers_postAge30, percActivate, totNumPpl, WHO_LEmale, SimuLEmale, WHO_LEage30_male, SimuLEage30_maleNonsmokers, SimuLEage30_maleSmokers, WHO_LEfemale, SimuLEfemale, WHO_LEage30_female, SimuLEage30_femaleNonsmokers, SimuLEage30_femaleSmokers';
    tablePrinter(rowLabelStr, meansAndCIs_1990, 'meansAndCIs_1990', validationGraphsFold);
    tablePrinter(rowLabelStr, meansAndCIs_2000, 'meansAndCIs_2000', validationGraphsFold);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if treatmentStatistics == 1
	%grab csv files
	cd(strcat(masterFold, baseCase))
	treatmentDemog = dlmread('treatmentDemog.csv', ',' ,1,0);
	diagnosedPeople = dlmread('diagnosedPpl.csv', ',' ,1,0);
    DOTSperfMeas = dlmread('DOTSperformanceMeas.csv', ',' ,1,0);
%     numDead_CatI = dlmread('NumberDead_CatIOnly.csv', ',' ,1,0);
%     numDead_CatII = dlmread('NumberDead_CatIIOnly.csv', ',' ,1,0);
  
	%frac of total pop registered this quarter with no prior treatment
	disp('frac of total pop registered this quarter with no prior treatment')
	fracRegisteredCATI_2010Q3 =  sum(diagnosedPeople(1734:1736, 6)) / diagnosedPeople(1735, 1)
	%frac of total pop registered this quarter with  prior treatment
	disp('frac of total pop registered this quarter with  prior treatment')
	fracRegisteredCATII_2010Q3 = sum(diagnosedPeople(1734:1736, 7)) / diagnosedPeople(1735, 1)
    %frac of total pop registered this quarter 
	disp('frac of total pop registered this quarter')
	fracRegisteredtot_2010Q3 = (sum(diagnosedPeople(1734:1736, 6)) + sum(diagnosedPeople(1734:1736, 7))) / diagnosedPeople(1735, 1)
	%catI vs catII registrations
	disp('catI vs catII registrations ')
	CatIvsCatIIregis2010Q3 = sum(diagnosedPeople(1734:1736, 6)) / sum(diagnosedPeople(1734:1736, 7))
    %initial num ppl in treatment (catI/II) + ppl diagnosed / total population
    disp('num ppl in treatment 2010Q3')
    numPplInTreatment2010Q3 = (treatmentDemog(175, 1) + treatmentDemog(175, 3) + treatmentDemog(175, 5) + treatmentDemog(175, 7) + sum(sum(diagnosedPeople(1734:1736, 6:7))) )/ diagnosedPeople(1735, 1)
    
	%treament demographics
	disp('--------------------------------')
	disp('CAT I proportion male ave2010')
	CatIPropMale2010 = mean(treatmentDemog(171:182, 1) ./ (treatmentDemog(171:182, 1) + treatmentDemog(171:182, 3)))
	
	disp('CAT II proportion male ave2010')
	CatIIPropMale2010 = mean(treatmentDemog(171:182, 5) ./ (treatmentDemog(171:182, 5) + treatmentDemog(171:182, 7)))
	
	disp('CAT I proportion old ave2010')
	CatIPropOld2010 = 1 - mean(  (treatmentDemog(171:182, 2) +  treatmentDemog(171:182, 4) ) ./ (treatmentDemog(171:182, 1) + treatmentDemog(171:182, 3)))
	
	disp('CAT II proportion old ave2010')
	CatIIPropOld2010 = 1 - mean(  (treatmentDemog(171:182, 6) +  treatmentDemog(171:182, 8) ) ./ (treatmentDemog(171:182, 5) + treatmentDemog(171:182, 7)))
    
    disp('catI demog balance')
    CatIPropYoungMale2010 = mean(treatmentDemog(171:182, 2) ./ (treatmentDemog(171:182, 1) + treatmentDemog(171:182, 3)))
	CatIPropYoungFemale2010 = mean(treatmentDemog(171:182, 4) ./ (treatmentDemog(171:182, 1) + treatmentDemog(171:182, 3)))	
    CatIPropOldMale2010 = CatIPropMale2010 - CatIPropYoungMale2010
	CatIPropOldFemale2010 = (1-CatIPropMale2010) - CatIPropYoungFemale2010
    disp('catII demog balance')
    CatIIPropYoungMale2010 = mean(treatmentDemog(171:182, 6) ./ (treatmentDemog(171:182, 5) + treatmentDemog(171:182, 7)))
	CatIIPropYoungFemale2010 = mean(treatmentDemog(171:182, 8) ./ (treatmentDemog(171:182, 5) + treatmentDemog(171:182, 7)))
    CatIIPropOldMale2010 = CatIIPropMale2010 - CatIIPropYoungMale2010
	CatIIPropOldFemale2010 = (1-CatIIPropOldMale2010) - CatIIPropYoungFemale2010
	
    
    disp('--------------------------------')
    %treatment quality 
%     disp('CAT I death rate')
%     CatIdeathRate_2010Q3 = sum(numDead_CatI(1737:1739,1)) / max(DOTSperfMeas(178:180, 1))
%     disp('CAT II death rate')
%     CatIIdeathRate_2010Q3 = sum(numDead_CatII(1737:1739,1)) / max(DOTSperfMeas(178:180, 5))
%     
    disp('CAT I default rate')
    CATIdefaults_2010Q3 = sum(DOTSperfMeas(178:180, 2))/max(DOTSperfMeas(178:180, 1))
    disp('CAT II default rate')
    CATIIdefaults_2010Q3 = sum(DOTSperfMeas(178:180, 6))/max(DOTSperfMeas(178:180, 5))
    
    disp('CAT I Failure rate')
    CATIfailures_2010Q3 = mean((DOTSperfMeas(178:180, 3) - DOTSperfMeas(178:180, 4)) ./ DOTSperfMeas(178:180, 3))
    disp('CAT II Failure rate')
    CATIIfailures_2010Q3 = mean((DOTSperfMeas(178:180, 7) - DOTSperfMeas(178:180, 8)) ./ DOTSperfMeas(178:180, 7))
    
    
    disp('--------------------------------')
    format('longG')
    disp('basic stats')
    [fracRegisteredCATI_2010Q3; fracRegisteredCATII_2010Q3; fracRegisteredtot_2010Q3; CatIvsCatIIregis2010Q3; CatIPropMale2010; CatIIPropMale2010; CatIPropOld2010; CatIIPropOld2010; CATIdefaults_2010Q3; CATIIdefaults_2010Q3 ]
    format('shortG')
    disp('catI balance')
    [CatIPropYoungMale2010; CatIPropYoungFemale2010; CatIPropOldMale2010; CatIPropOldFemale2010]
    disp('catII balance')
    [CatIIPropYoungMale2010; CatIIPropYoungFemale2010; CatIIPropOldMale2010; CatIIPropOldFemale2010]
    
end
