function resultsPlotter(recordUnit, folderName, plotResolution)

%this plotter makes the graphs for comparing sex-age and smoking-age

%prevalences for lat DS, lat MDR, act DS, act MDR

%recordUnit is in months

%folderName is the folder that has the outcome graphs

%plot resolution can be '-r70'



%plot to make results graphs

%ageGrp1 = 0_20

%ageGrp2 = 21_60

%ageGrp3 = 61_100



%read in .csv files and get data

currentDirectory = pwd;

cd(folderName);  %go to the output directory



female = dlmread('HealthOutcomes_femaleOnly.csv', ',' ,1,6);

male = dlmread('HealthOutcomes_maleOnly.csv', ',' ,1,6);

male_ageGrp1 = dlmread('HealthOutcomes_male_age0_20Only.csv', ',' ,1,6);

male_ageGrp2 = dlmread('HealthOutcomes_male_age21_60Only.csv', ',' ,1,6);

male_ageGrp3 = dlmread('HealthOutcomes_male_age61_100Only.csv', ',' ,1,6);

female_ageGrp1 = dlmread('HealthOutcomes_female_age0_20Only.csv', ',' ,1,6);

female_ageGrp2 = dlmread('HealthOutcomes_female_age21_60Only.csv', ',' ,1,6);

female_ageGrp3 = dlmread('HealthOutcomes_female_age61_100Only.csv', ',' ,1,6);



smoking_ageGrp1 = dlmread('HealthOutcomes_smoking_age0_20Only.csv', ',' ,1,6);

smoking_ageGrp2 = dlmread('HealthOutcomes_smoking_age21_60Only.csv', ',' ,1,6);

smoking_ageGrp3 = dlmread('HealthOutcomes_smoking_age61_100Only.csv', ',' ,1,6);

nonSmoking_ageGrp1 = dlmread('HealthOutcomes_nonSmoking_age0_20Only.csv', ',' ,1,6);

nonSmoking_ageGrp2 = dlmread('HealthOutcomes_nonSmoking_age21_60Only.csv', ',' ,1,6);

nonSmoking_ageGrp3 = dlmread('HealthOutcomes_nonSmoking_age61_100Only.csv', ',' ,1,6);



noSmokeDeathAgeFreq = dlmread('nosmokeDeathAgesMat.csv', ',' ,1,0); 

smokerDeathAgeFreq = dlmread('smokeDeathAgesMat.csv', ',' ,1,0);

noSmokerDSTBDeathAgeFreq  = dlmread('nosmokeDSTBdeathAgesMat.csv', ',' ,1,0); 

smokerDSTBDeathAgeFreq = dlmread('smokeDSTBdeathAgesMat.csv', ',' ,1,0);

noSmokerMDRTBDeathAgeFreq = dlmread('nosmokeMDRTBdeathAgesMat.csv', ',' ,1,0);

smokerMDRTBDeathAgeFreq = dlmread('smokeMDRTBdeathAgesMat.csv', ',' ,1,0);



cd(currentDirectory); %get back out to the original directory



%%%%%%%%%%%%%% data captured %%%%%%%%%%%%%%%%





%%%%%%%%%%%% smoking ratio graph %%%%%%%%%%%

%data is sparse so bin across years (5 year bins, 5 year ages)



matLength = size(noSmokeDeathAgeFreq,1);

matWidth = size(noSmokeDeathAgeFreq,2);



noSmokeDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

noSmokeDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);

smokerDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

smokerDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);

noSmokerDSTBDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

noSmokerDSTBDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);

smokerDSTBDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

smokerDSTBDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);

noSmokerMDRTBDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

noSmokerMDRTBDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);

smokerMDRTBDeathAgeFreqShort = zeros( floor((matLength-5)/5)+1  , size(noSmokeDeathAgeFreq,2));

smokerMDRTBDeathAgeFreqSmall = zeros(floor((matLength-5)/5)+1 , floor((matWidth-5)/5)+1);



for x = 0:floor((matLength-5)/5)

    noSmokeDeathAgeFreqShort(x+1,:) = sum(noSmokeDeathAgeFreq(5*x+1:5*x+5,:));

    smokerDeathAgeFreqShort(x+1,:) = sum(smokerDeathAgeFreq(5*x+1:5*x+5,:));

    noSmokerDSTBDeathAgeFreqShort(x+1,:) = sum(noSmokerDSTBDeathAgeFreq(5*x+1:5*x+5,:));

    smokerDSTBDeathAgeFreqShort(x+1,:) = sum(smokerDSTBDeathAgeFreq(5*x+1:5*x+5,:));

    noSmokerMDRTBDeathAgeFreqShort(x+1,:) = sum(noSmokerMDRTBDeathAgeFreq(5*x+1:5*x+5,:));

    smokerMDRTBDeathAgeFreqShort(x+1,:) = sum(smokerMDRTBDeathAgeFreq(5*x+1:5*x+5,:));

end

for x = 0:floor((matWidth-5)/5)

    noSmokeDeathAgeFreqSmall(:, x+1) = sum( (noSmokeDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

    smokerDeathAgeFreqSmall(:, x+1) = sum( (smokerDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

    noSmokerDSTBDeathAgeFreqSmall(:, x+1) = sum( (noSmokerDSTBDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

    smokerDSTBDeathAgeFreqSmall(:, x+1) = sum( (smokerDSTBDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

    noSmokerMDRTBDeathAgeFreqSmall(:, x+1) = sum( (noSmokerMDRTBDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

    smokerMDRTBDeathAgeFreqSmall(:, x+1) = sum( (smokerMDRTBDeathAgeFreqShort(:,5*x+1:5*x+5)') )';

end



noTBsmokeRatio = smokerDeathAgeFreqSmall./noSmokeDeathAgeFreqSmall;

DsTBsmokeRatio = smokerDSTBDeathAgeFreqSmall./noSmokerDSTBDeathAgeFreqSmall;

MdrTBsmokeRatio = smokerMDRTBDeathAgeFreqSmall./noSmokerMDRTBDeathAgeFreqSmall;



noTBsmokeRatio(isnan(noTBsmokeRatio)) = 0;

DsTBsmokeRatio(isnan(DsTBsmokeRatio)) = 0;

MdrTBsmokeRatio(isnan(MdrTBsmokeRatio)) = 0;



ageVec = [1:5:100];

plot(ageVec, noTBsmokeRatio(38,:), 'b');  %38 5-yr periods after burn in, or roughly 2012

hold on;

plot(ageVec, DsTBsmokeRatio(38,:), 'r');

hold on;

plot(ageVec, MdrTBsmokeRatio(38,:), 'g');

leg = legend('no TB', 'DS TB','MDR TB');

set(leg,'Location','NorthEastOutside');

xlabel('age');

ylabel('ratio of smoker deaths to nonsmokers deaths');

title('ratio of smoking to nonsmoking deaths, by age');

print('-dpng',plotResolution,[folderName '/' 'smokingDeathAgesRatioGraph']);



%%%%%%%%%prop graphs%%%%%%%%%%%



%make x-axis labels

numEntries = size(female, 1);

startTreatmentIndex = (130*12)/recordUnit;  %130 year burn in period, record unit is in units of month

%endTreatmentIndex = startTreatmentIndex+((12*40)/recordUnit); %40 is 15 years from 1996-2011 + 25 years projection
endTreatmentIndex = startTreatmentIndex+((12*29)/recordUnit); %29 is to 2025

%xCoord_postBurnIn = [1996:(recordUnit/12):1996+40];
xCoord_postBurnIn = [1996:(recordUnit/12):1996+29];


titlePrefixArray = {'Proportion of Demographic with Latent DS: ','Proportion of Demographic with Latent MDR: ','Proportion of Demographic with Active DS: ','Proportion of Demographic with Active MDR: '};

pngtTitlePrefixArray = {'PropLatDS_','PropLatMDR_','PropActDS_','PropActMDR_'};



%PLOT

%for loop here to loop over latent DS, latent MDR, active DS, active MDR (cols 6, 7, 8, and 9 of the tables, cols 1 to 4 of the dlm matracies)

for outcomeIndex = 1:4

    

    titlePrefix = titlePrefixArray{outcomeIndex};

    pngTitlePrefix = pngtTitlePrefixArray{outcomeIndex};

    

    %sex and age

    plot(xCoord_postBurnIn, male_ageGrp1(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'b' );

    hold on;

    plot(xCoord_postBurnIn, female_ageGrp1(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'r' );

    hold on;

    plot(xCoord_postBurnIn, male_ageGrp2(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'b-.' );

    hold on;

    plot(xCoord_postBurnIn, female_ageGrp2(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'r-.' );

    hold on;

    plot(xCoord_postBurnIn, male_ageGrp3(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'b:' );

    hold on;

    plot(xCoord_postBurnIn, female_ageGrp3(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'r:' );



    leg = legend('male 0-20 yrs', 'female 0-20 yrs','male 21-60 yrs','female 21-60 yrs', 'male 61-100 yrs', 'female 61-100 yrs');

    set(leg,'Location','NorthEastOutside');

    titleSt = strcat(titlePrefix, 'Sex and Age');

    title(titleSt);

    xlabel('Year')

    ylabel('Proportion of Demographic Group')

    pngTitle = strcat(pngTitlePrefix, 'Sex and Age'); 

    axis([1996 (1996+40) 0 1])

    axis 'auto y'

    print('-dpng',plotResolution,[folderName '/' pngTitle]);

    close all



    %smoking and age

    plot(xCoord_postBurnIn, smoking_ageGrp1(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'k' );

    hold on;

    plot(xCoord_postBurnIn, nonSmoking_ageGrp1(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'m' );

    hold on;

    plot(xCoord_postBurnIn, smoking_ageGrp2(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'k-.' );

    hold on;

    plot(xCoord_postBurnIn, nonSmoking_ageGrp2(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'m-.' );

    hold on;

    plot(xCoord_postBurnIn, smoking_ageGrp3(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'k:' );

    hold on;

    plot(xCoord_postBurnIn, nonSmoking_ageGrp3(startTreatmentIndex:endTreatmentIndex, outcomeIndex),'m:' );



    leg = legend('smokers 0-20 yrs', 'nonsmokers 0-20 yrs', 'smokers 21-60 yrs', 'nonsmokers 21-60 yrs', 'smokers 61-100 yrs', 'nonsmokers 61-100 yrs');

    set(leg,'Location','NorthEastOutside');

    titleSt = strcat(titlePrefix, 'Smoking Status and Age');

    title(titleSt);

    xlabel('Year');

    ylabel('Proportion of Demographic Group');

    pngTitle = strcat(pngTitlePrefix, 'Smoking Status and Age'); 

    axis([1996 (1996+40) 0 1])

    axis 'auto y'

    print('-dpng',plotResolution,[folderName '/' pngTitle]);

    close all



end %end for loop



