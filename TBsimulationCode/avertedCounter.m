function avertedCounter(folderNameTreatment, treatmentDescrip, folderNameBaseline, baselineDescrip, plotResolution, folderName)
%this function makes deaths averted and disease cases averted tables/graphs
%recordUnit is in months
%folderNameTreatment is the folder that has the treatment outcome graphs
%treatmentDescrip is a string to describe the treatment
%folderNameTreatment is the folder that has the baseline outcome graphs (i.e., no treatment folder)
%baselineDescrip is a string to describe the baseline
%plot resolution can be '-r70'

%plot to make results graphs
%ageGrp1 = 0_20
%ageGrp2 = 21_60
%ageGrp3 = 61_100

%read in .csv files and get data
currentDirectory = pwd;
cd(folderNameTreatment);  %go to the output directory

treatmentDiseaseCases = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);   
treatmentDeaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
treatmentQalys = dlmread('totalQalys.csv', ',' ,1,0);

cd(folderNameBaseline);  %go to the output directory

baselineDiseaseCases = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);  
baselineDeaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
baselineQalys = dlmread('totalQalys.csv', ',' ,1,0);

cd(currentDirectory); %get back out to the original directory
mkdir(folderName);
cd(folderName);

%PLOT

latDSinfectionsAverted =  treatmentDiseaseCases(193:481,1) - baselineDiseaseCases(193:481,1);  %193:481 is 2012 to 2036 inclusive, ie 25 yrs projection
latMDRinfectionsAverted =  treatmentDiseaseCases(193:481,2) - baselineDiseaseCases(193:481,2);
actDScasesAverted = treatmentDiseaseCases(193:481,3) - baselineDiseaseCases(193:481,3);
actMDRcasesAverted = treatmentDiseaseCases(193:481,4) - baselineDiseaseCases(193:481,4);

recordUnit = 1;
yrs = [2012:(recordUnit/12):(1996+40)];  %40 years after 1996 is also equal to 25 years projection from jan 2012

subplot(2,2,1), plot(yrs, latDSinfectionsAverted);
text(2015,(max(latDSinfectionsAverted) + 500),{'Latent DS Infections Averted'},'VerticalAlignment','Bottom');
xlabel('Year');
ylabel('Number Averted');
subplot(2,2,2), plot(yrs, actDScasesAverted);
text(2015,(max(actDScasesAverted) + 30),{'Active DS Cases Averted'},'VerticalAlignment','Bottom');
xlabel('Year');
ylabel('Number Averted');
subplot(2,2,3), plot(yrs, latMDRinfectionsAverted);
text(2015,(max(latMDRinfectionsAverted) + 50),{'Latent MDR Infections Averted'},'VerticalAlignment','Bottom');
xlabel('Year');
ylabel('Number Averted');
subplot(2,2,4), plot(yrs, actMDRcasesAverted);
text(2015,(max(actMDRcasesAverted) + 3),{'Active MDR Cases Averted'},'VerticalAlignment','Bottom');
titleStr = strcat(treatmentDescrip, ' compared to ', baselineDescrip);
xlabel('Year');
ylabel('Number Averted');

title(titleStr);
print('-dpng',plotResolution,[currentDirectory '/' 'CasesAverted']);
close all

totLatDSinfectionsAverted = sum(latDSinfectionsAverted);
totLatMDRinfectionsAverted =  sum(latMDRinfectionsAverted);
totActDScasesAverted = sum(actDScasesAverted);
totActMDRcasesAverted = sum(actMDRcasesAverted);
totQalys = treatmentQalys-baselineQalys;

totalAverted = [totLatDSinfectionsAverted, totLatMDRinfectionsAverted, totActDScasesAverted, totActMDRcasesAverted, totQalys];  
tableStr = strcat(treatmentDescrip, '_vs_', baselineDescrip);
tablePrinter('totLatDSinfectionsAverted, totLatMDRinfectionsAverted, totActDScasesAverted, totActMDRcasesAverted, totQalys', totalAverted, tableStr, folderName);


