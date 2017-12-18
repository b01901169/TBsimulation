%calibration Target maker
clear
clear all
close all

masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\May28\';
baseFolder = strcat(masterFolder,'base\');
currentDirectory = pwd;

%%use this if many calibration runs you want to do this on: all in base
% foldersStruc = dir(masterFolder);
% is_dir = [foldersStruc(:).isdir]';
% foldnameArray_all = {foldersStruc(:).name}';
% foldNameArray = foldnameArray_all(is_dir);
% foldNameArray = foldNameArray(3:end); % Exclude . and ..

%use this if just one base to do it on (latest run)
cd(masterFolder);
folderNameStartsWith = 'p01';  %this is the first letter of the runs
mustHaveCSV = 'diagnosisTime_CatI.csv';
foldNameArray{1} = findMostRecentLegitFolder(baseFolder, folderNameStartsWith, mustHaveCSV);
foldNameArray{1} = strrep(foldNameArray{1}, masterFolder, '');
cd(currentDirectory)
    

%manual list
% foldNameArray = {
%     's01_2014-04-03_13-24-14'
%     's01_2014-04-03_14-40-09'
%     's01_2014-04-03_15-51-11'
%     's01_2014-04-03_17-01-49'
%     };

for i = 1:size(foldNameArray,1)
%for i = 1:1
    
    folderName = strcat(masterFolder, foldNameArray{i});
%    folderName = masterFolder;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %grab data
    
    cd(folderName);  %go to the output directory
    
    %TBsimulation.txt
    string = textread('TB simulation.txt', '%s');
    
    %time till diag
    diagTime = dlmread('diagnosisTime_CatI.csv', ',' ,1,0);
    
    %healthOutcomes
    healthOut = dlmread('HealthOutcomes.csv', ',' ,1,0);
    
    %validationGraphMaker stuff
	treatmentDemog = dlmread('treatmentDemog.csv', ',' ,1,0);
	diagnosedPeople = dlmread('diagnosedPpl.csv', ',' ,1,0);
    DOTSperfMeas = dlmread('DOTSperformanceMeas.csv', ',' ,1,0);
%    numDead_CatI = dlmread('NumberDead_CatIOnly.csv', ',' ,1,0);
%    numDead_CatII = dlmread('NumberDead_CatIIOnly.csv', ',' ,1,0);
    
    %age structure of prev
    %ageStrucPrev{1} = dlmread('hasTBin_1525.csv', ',' ,1,0);  %1993
    ageStrucPrev{1} = dlmread('hasTBin_1585.csv', ',' ,1,0);  %1998
    ageStrucPrev{2} = dlmread('hasTBin_1645.csv', ',' ,1,0);  %2003
    ageStrucPrev{3} = dlmread('hasTBin_1669.csv', ',' ,1,0);  %2005
    
    %prev graph
	simOut = dlmread('monthlyActOutcomes.csv', ',' ,1,0);
    cd(currentDirectory); %get back out to the original directory
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %latToAct_cal
    c_latToAct_cal = str2num(string{15});
    
    %FOI_cal
    c_FOI_cal = str2num(string{19});
    
    %aveUptake_cal
    c_aveUptake_cal =str2num(string{23});
    
    %cat2uptake_cal
    c_cat2uptake_cal =str2num(string{27});
    
    %ave time till uptake
    c_diagTime = diagTime(1737,1);
    
    %health outcomes (latent and act at 1996)
    c_healthOut = healthOut(25,[7 9]);
    

	%frac of total pop registered this quarter with no prior treatment
	fracRegisteredCATI_2010Q3 =  sum(diagnosedPeople(1734:1736, 6)) / diagnosedPeople(1735, 1)

	%frac of total pop registered this quarter with  prior treatment
	fracRegisteredCATII_2010Q3 = sum(diagnosedPeople(1734:1736, 7)) / diagnosedPeople(1735, 1)

    %frac of total pop registered this quarter 
	fracRegisteredtot_2010Q3 = (sum(diagnosedPeople(1734:1736, 6)) + sum(diagnosedPeople(1734:1736, 7))) / diagnosedPeople(1735, 1)

	%catI vs catII registrations
	CatIvsCatIIregis2010Q3 = sum(diagnosedPeople(1734:1736, 6)) / sum(diagnosedPeople(1734:1736, 7))

    %initial num ppl in treatment (catI/II) + ppl diagnosed / total population
    numPplInTreatment2010Q3 = (treatmentDemog(175, 1) + treatmentDemog(175, 3) + treatmentDemog(175, 5) + treatmentDemog(175, 7) + sum(sum(diagnosedPeople(1734:1736, 6:7))) )/ diagnosedPeople(1735, 1)
    
    valOut(i,:) = [fracRegisteredCATI_2010Q3, fracRegisteredCATII_2010Q3, fracRegisteredtot_2010Q3, CatIvsCatIIregis2010Q3];
    summaryMat(i,:) = [c_latToAct_cal, c_FOI_cal, c_aveUptake_cal, c_cat2uptake_cal, c_diagTime, c_healthOut];
    
    %make prev age structure graphs
    ages = [0	5	10	15	20	25	30	35	40	45	50	55	60	65	70	75	80	85	90	95];
    indiaData = 100000*[
        0.000885419	0.001385228	0.001258286	0.001434504	0.002603419	0.003761941	0.00555499	0.006804255	0.010593589	0.008126246	0.009727112	0.011923014	0.01566252	0.015617919	0.016139674	0.013221792	0.013140364	0.011344524	0.012633881	0.006784333
        0.001501273	0.001531214	0.001455726	0.002549121	0.003697059	0.005693879	0.007811364	0.009402176	0.008277181	0.010608352	0.009335829	0.012405081	0.01182184	0.013570467	0.015404721	0.014522619	0.013048252	0.012983076	0.003347435	0.0249686
        0.010106094	0.001261491	0.001601428	0.001588536	0.002062728	0.002241193	0.003290454	0.004530522	0.005223636	0.005470236	0.006642283	0.007415469	0.007393865	0.006237989	0.007526925	0.007147744	0.006105058	0.006387843	0.007832773	0.003363833
        0.000804808	0.001233151	0.001081742	0.002420975	0.002657252	0.004153483	0.005391537	0.006259815	0.007221397	0.00845402	0.010197851	0.008638978	0.010048536	0.011496754	0.010979185	0.007245091	0.007496174	0.004508902	0.002765955	0.006865667
        ];
    indiaDataLCI = 100000*[
        0.000564831	0.001007675	0.000845952	0.00093709	0.001884396	0.002907595	0.004329588	0.005436281	0.008442985	0.00627512	0.007544534	0.009349361	0.012512291	0.011730899	0.0119059	0.007795961	0.007940223	0.003811574	0.003933493	0.001143467
        0.001010674	0.001100023	0.001039515	0.001884944	0.002839335	0.004609168	0.006467801	0.00790952	0.006646203	0.008474209	0.007250663	0.009536031	0.009196346	0.010322229	0.011456772	0.009042127	0.007235848	0.005199902	0.000473309	0.0067356
        0.002122819	0.000908963	0.001309402	0.001242699	0.00165252	0.001844381	0.002760593	0.003632985	0.004233842	0.004495699	0.005302844	0.00573152	0.005981113	0.005010283	0.005741436	0.004543289	0.003851942	0.002221384	0.004602893	0.000942
        0.000419211	0.000782717	0.000670554	0.001668544	0.001889094	0.003170265	0.004132469	0.004906966	0.005586225	0.006580491	0.007802262	0.006440876	0.007572606	0.008162713	0.00761641	0.00364709	0.003971685	0.001588077	0.000729716	0.0011611
        ];
    
    indiaDataUCI = 100000*[
        0.001390832	0.001904312	0.00187734	0.002196598	0.003596382	0.004867075	0.007124879	0.00851363	0.013294399	0.010591727	0.012571944	0.015208158	0.019619489	0.020843446	0.021958049	0.022835361	0.022415829	0.033266792	0.042131291	0.039282433
        0.002231509	0.002131435	0.002039953	0.003447537	0.004812714	0.007032164	0.00943545	0.011177686	0.010312717	0.013280662	0.012050195	0.016123214	0.01521887	0.017888136	0.020730363	0.023251284	0.023477182	0.037844902	0.023100967	0.088155033
        0.055653297	0.001766028	0.001958473	0.00203041	0.002574925	0.002723163	0.003921618	0.005649594	0.006443699	0.006654975	0.008319487	0.009588911	0.009241577	0.007764778	0.009870674	0.011695825	0.009770294	0.018440139	0.013327118	0.014102633
        0.001553351	0.001955946	0.001752847	0.003513239	0.00374005	0.005442003	0.007031554	0.00799386	0.009355376	0.010927774	0.01334121	0.011634638	0.013363353	0.016255372	0.01601471	0.014696432	0.01439295	0.012833875	0.010473935	0.039483433
        ];
    
    distToUpperLimit = indiaDataUCI - indiaData;
    distToLowerLimit = indiaData - indiaDataLCI;
    
    %yearMat = [1993, 1998, 2003, 2005];
    yearMat = [1998, 2003, 2005];
    for i = 1:3
        simAgePrev = sum(ageStrucPrev{i}(:,3:4),2) ./ sum(ageStrucPrev{i}(:,1:2),2) * 100000;
        p = plot(ages, simAgePrev);
        set(p,'Color','red','LineWidth',1);
        xlabel('Age');
        ylabel('Active TB Prevalence (out of 100000 People)');
        hold on
        errorbar(ages, indiaData(i,:),distToLowerLimit(i,:),distToUpperLimit(i,:));
        graphName = strcat('TBprevNFHSComparison_', num2str(yearMat(i)));
        print('-dpng','-r70',[folderName '/' graphName]);
        ylim([0 9000])
        close all
    end
    
    %make a more resolved prev graph
	%PREVALENCE
    %who prevalence data
	dataYear = [1990    ; 1991  ; 1992  ; 1993  ; 1994  ;  1995 ; 1996  ; 1997  ; 1998  ; 1999  ; 2000  ; 2001  ; 2002  ; 2003  ; 2004  ; 2005  ; 2006  ; 2007  ; 2008  ; 2009  ; 2010  ];
	actTBPrev = [459 ; 460 ; 460    ; 461   ; 462   ; 462   ; 463   ; 464   ; 464   ; 465   ; 466   ; 466   ; 436   ; 409   ; 383   ; 358   ; 335   ; 314   ; 294   ; 275   ; 256];
	distToUpperLimit = [56  ; 56    ; 56    ; 56    ; 56    ; 57    ; 56    ; 56    ; 57    ; 56    ; 56    ; 57    ; 62    ; 66    ; 71    ; 78    ; 84    ; 91    ; 99    ; 107   ; 117];
	distToLowerLimit = [52  ; 53    ; 52    ; 53    ; 53    ; 53    ; 53    ; 53    ; 53    ; 53    ; 54    ; 53    ; 57    ; 62    ; 66    ; 70    ; 74    ; 80    ; 85    ; 90    ; 95];
	%simulation data
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
    print('-dpng','-r70',[folderName '/' 'TBprevWHOcomparisonFine']);
end

disp('targets are:')
disp('Proportion of population with latent TB, 1999: 0.44 from Dye')
disp('Proportion of population with active TB, 1999: 0.005 from Dye')
disp('fracRegisteredCATI_2010Q3 = 0.000263')
disp('fracRegisteredCATII_2010Q3 = 0.000064')
disp('fracRegisteredtot_2010Q3 = 0.000313')
disp('CatIvsCatIIregis2010Q3 = 4.07')
disp('WHO prevalence over time')
disp('WHO incidence over time')
disp('age stratified prevalence')
disp('MDR in 2008')

disp('folders are')
foldNameArray

disp('latToAct, FOI, aveUptake, cat2Uptake, diagnosis time in months, lat DS in 1996, act DS in 1996')
format
summaryMat

disp('frac of total pop registered this quarter with no prior treatment, with prior trt, total, and catI vs catII registrations')
format('longG')
valOut