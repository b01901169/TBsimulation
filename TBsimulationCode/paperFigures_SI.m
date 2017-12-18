%%SI graphs for PNAS


%FOLDERS
masterFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Aug10\';
baseCase = 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta\2012-08-10_21-40-28';
plotResolution = '-r100';

outputFolder = strcat(masterFold, 'paperFigures\SI_figures');
plotFixerFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs'; %this is where plotfixer lives;
currentDirectory = masterFold;

cd(masterFold);
TBsimParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Age structure figure
cd(strcat(masterFold,baseCase));
ageStrucMat = dlmread('AgeStructure.csv', ',' ,1,0);
cd(masterFold);
toPlot = ageStrucMat( 1:35, 6:10);
area([1871:5:2041], toPlot);
%hold on;
%plot([2013 2013], [0 1], 'w--');
%hold;
%plot([2013 2013], [0 1], 'w:');
ylim([0 1]);
xlabel('Year');
ylabel('Proportion of Total');
legend('20 and under', '21 to 40', '41 to 60', '61 to 80', '81 and above', 'Location','NorthWest');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'ageStructure.png']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoking mortality graph
ageGrps = TBparams.baseMortAgeBrac(:,1);

%male %female in 2000
overallMort =[...
0.005831265	,	0.005960499	;...
0.000454064	,	0.000697257	;...
0.000155821	,	0.00019748	;...
0.000117493	,	0.000121659	;...
0.000154155	,	0.000200813	;...
0.000213311	,	0.000259133	;...
0.000269964	,	0.00027163	;...
0.000323281	,	0.000259133	;...
0.000427409	,	0.00028246	;...
0.000538188	,	0.000344107	;...
0.000777198	,	0.000486548	;...
0.001147674	,	0.000770536	;...
0.001827495	,	0.001273355	;...
0.00255257	,	0.00183914	;...
0.004122313	,	0.003183256	;...
0.005508937	,	0.004633399	;...
0.008489591	,	0.007120363	;...
0.010384868	,	0.009100834	;...
0.013383796	,	0.012170335	;...
0.018170391	,	0.017015244	;...
0.025971778	,	0.024864003	;...
0.039035201	,	0.037926249	;...
];

%using 2000 parameters
subplot(2,1,1); %males
plot(ageGrps(1:end-6), mortMaleSmokingOverTime(1:end-6,2), 'r');
hold on;
plot(ageGrps(1:end-6), mortMaleNonsmokingOverTime(1:end-6,2), 'b');
hold on;
plot(ageGrps(1:end-6), overallMort(1:end-6,1), 'k');
title('Male');
xlabel('Age');
ylabel('Monthly Probability of Death');
%legend('Overall Mortality', 'Smokers', 'Nonsmokers','Location','NorthWest');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);

subplot(2,1,2); %females
plot(ageGrps(1:end-6), mortFemaleSmokingOverTime(1:end-6,2), 'r');
hold on;
plot(ageGrps(1:end-6), mortFemaleNonsmokingOverTime(1:end-6,2), 'b');
hold on;
plot(ageGrps(1:end-6), overallMort(1:end-6,2), 'k');
title('Female');
xlabel('Age');
ylabel('Monthly Probability of Death');
legend('Smokers', 'Nonsmokers','Overall Mortality', 'Location','NorthWest');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'smokingMort.png']);
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smoking start and stop rates (male, by smoking duration)
urbFrac = TBparams.urbanFrac(1);
rurFrac = 1- TBparams.urbanFrac(1);
zero = urbFrac*TBparams.oneMonQuittingUrbanProb*(TBparams.smoker_quit_ratio^0) +  rurFrac*TBparams.oneMonQuittingRuralProb*(TBparams.smoker_quit_ratio^0)  ;
one = urbFrac*TBparams.oneMonQuittingUrbanProb*(TBparams.smoker_quit_ratio^1)  + rurFrac*TBparams.oneMonQuittingRuralProb*(TBparams.smoker_quit_ratio^1);
two = urbFrac*TBparams.oneMonQuittingUrbanProb*(TBparams.smoker_quit_ratio^2) + rurFrac*TBparams.oneMonQuittingRuralProb*(TBparams.smoker_quit_ratio^2) ;
three = urbFrac*TBparams.oneMonQuittingUrbanProb*(TBparams.smoker_quit_ratio^3) +  rurFrac*TBparams.oneMonQuittingRuralProb*(TBparams.smoker_quit_ratio^3) ;
four = urbFrac*TBparams.oneMonQuittingUrbanProb*(TBparams.smoker_quit_ratio^4)  + rurFrac*TBparams.oneMonQuittingRuralProb*(TBparams.smoker_quit_ratio^4);

plot(TBparams.smokingChurnAgeBrac(1:65,1), (zero(1:65,1)), 'Color', [0 0 1]);
hold on;
plot(TBparams.smokingChurnAgeBrac(1:65,1), (one(1:65,1)),'Color', [0 0.2 0.9]);
hold on;
plot(TBparams.smokingChurnAgeBrac(1:65,1), (two(1:65,1)),'Color', [0 0.4 0.8]);
hold on;
plot(TBparams.smokingChurnAgeBrac(1:65,1), (three(1:65,1)),'Color', [0 0.6 0.7]);
hold on;
plot(TBparams.smokingChurnAgeBrac(1:65,1), (four(1:65,1)),'Color', [0 0.8 0.6]);
xlabel('Age');
ylabel('Probability of Quitting Smoking');
legend('<1 year of smoking','1-2 years of smoking','2-3 years of smoking','3-4 years of smoking','4 or more years of smoking', 'Location','NorthEast');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'smokingChurn.png']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOTS coverage over time
plot(empiricalYear, empiricalCoverage, 'o');
hold on;
plot(coverageNeededAt, TBparams.DotsCoverageSequence);
xlabel('Year');
ylabel('Proportion of India Covered by DOTS');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'DOTScoverage.png']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% treatment uptake levels by age, sex and treatment category  
urban = TBparams.prTestedGivenTB*TBparams.SSsensit * TBparams.uptakeUrbanKnowledge;
rural = TBparams.prTestedGivenTB*TBparams.SSsensit * TBparams.uptakeRuralKnowledge;
cat2 =  TBparams.prTestedGivenTB*TBparams.SSsensit* TBparams.catIIuptakeKnowledge .* TBparams.priorTreatBoostFactor;

getTreatment =(repmat(TBparams.urbanFrac,size(TBparams.uptakeUrbanKnowledge,1),1).*urban) + (( repmat(1-TBparams.urbanFrac,size(TBparams.uptakeUrbanKnowledge,1),1) ) .* rural );

subplot(2,1,1);
plot(TBparams.uptakeAgeBracs(:,1), getTreatment);
xlabel('Age');
ylabel({'Probability of Entering Cat. I/III','Treatment If No Treatment History'});
legend('Male', 'Female', 'Location','NorthWest');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);

subplot(2,1,2);
plot(TBparams.uptakeAgeBracs(:,1), cat2);
xlabel('Age');
ylabel({'Probability of Entering Cat. II','Treatment If Have Prior Treatment'});
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'getTreatment.png']);
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DOTS+ ramp up
%ramp up over 8 years (this is a 96 x 1 vector)  This was fitted to an exp
%using the rntcpt report dots plus data
plot([2007:1/12:2015-(1/12)] , TBparams.DotsPlusRampUpSequence);
xlabel('Years');
ylabel('DOTS-Plus Coverage Over Time');
cd(plotFixerFold);
plotFixerPNAS;
cd(currentDirectory);
print('-dpng',plotResolution,[outputFolder '/' 'DotsPlusRampUp.png']);
close all






