%making graphs for RIPS ppt

plotResolution = '-r250';

startYear = 26; %this is starting in 1996
%startYear = 29; %this is starting in 2011

%%%%%%%NOTES:
% % Note that the no cat4 / ave cat4 / perfect cat4 graphs can have the "perfect dots" above the non-perfect dots for the no cat4 case
% % because "perfect dots" still has failures (people who graduate but have latent DS/MDR and also MDR people still fail).  The more
% % people you push to graduate the higher the number of people who graduate with latent (since it's just a fixed frac of those who
% % graduate) and therefore higher MDR.  To fix this: can make "perfect DOTS" mean that not only are the patients perfectly behaved
% % but also the SS test is perfect at finding graduates who would leave with latent and would catch all the MDR graduates.  This has
% % not been done.
%%%%%%%%%%%


%initialize folders
plotFixerFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs'; %this is where plotfixer lives;

%grab the data
masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\June27\';
pptGraphsFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\June27\validationGraphs'; %this is where the graps should be outputted to

%need to also add the subfolder names, once they are generated
aveTreat_fullCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta\2012-06-27_15-02-25');
% aveTreat_superfastDST = strcat(masterFolder, 'inf0p0020_lat2p55_base_Post2011_aveTreat_fullCatIV_empUptake_superfastDST\2012-02-03_04-01-44');  %replace with super fast
noEvolve = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve\2012-07-02_18-53-03');

% %aveTreat_noCat4_ultraUptake = strcat(masterFolder, 'inf0p0020_lat2p55_base_Post2011_aveTreat_fullCatIV_ultraUptake\2012-01-30_02-13-11');
aveTreat_noCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_aveTreat_noCatIV_empUpta\2012-06-29_07-30-41');
aveTreat_perfectCat4_empUptake = strcat(masterFolder,'inf0p0023_lat2p16_base_Post2011_aveTreat_perfectCatIV_empUpta\2012-06-29_09-21-55');

bestTreat_fullCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_bestTreat_fullCatIV_empUpta\2012-06-29_00-06-57');
bestTreat_noCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_bestTreat_noCatIV_empUpta\2012-06-29_11-16-48');
bestTreat_perfectCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_bestTreat_perfectCatIV_empUpta\2012-06-29_13-11-38');

perfectTreat_fullCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_perfectTreat_fullCatIV_empUpta\2012-06-29_15-07-31');
perfectTreat_noCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_perfectTreat_noCatIV_empUpta\2012-06-29_17-09-30');
perfectTreat_perfectCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_perfectTreat_perfectCatIV_empUp\2012-06-29_19-08-06');

worstTreat_fullCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_worstTreat_fullCatIV_empUpta\2012-06-29_21-05-27');
worstTreat_noCat4_empUptake = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_worstTreat_noCatIV_empUpta\2012-06-29_22-58-23');

noTreat = strcat(masterFolder, 'inf0p0023_lat2p16_noTreat_fullCatIV_empUpta\2012-06-30_00-51-46');
%
% %The master folders for these have to be put in manually as they are different than the rest
aveTreat_noEvolve_2019start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve2019\2012-07-02_20-47-52');
aveTreat_noEvolve_2027start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_noEvolve2027\2012-07-02_22-42-31');

bestTreat_2019start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_bestTreat2019\2012-06-30_21-20-18');
bestTreat_2027start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_bestTreat2027\2012-06-30_23-12-01');

%the half DOTS+ noEvolve scenarios
halfDotsPlus = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta\2012-07-05_19-22-25');
halfDotsPlusNoEVOLVE_2013start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta_noEvolve2013\2012-07-06_01-00-15');
halfDotsPlusNoEVOLVE_2019start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta_noEvolve2019\2012-07-05_23-28-34');
halfDotsPlusNoEVOLVE_2027start = strcat(masterFolder,'inf0p0023_lat2p16_aveTreat_halfCatIV_empUpta_noEvolve2027\2012-07-05_19-41-23');

% %the DST scenarios
aveTreat_allGeneXscenarios = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_allGeneXscenarios\2012-06-30_04-38-55');
aveTreat_GeneXinitialDST = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_GeneXinitialDST\2012-06-30_08-22-10');
aveTreat_GeneXinsteadOfLJ = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_GeneXinsteadOfLJ\2012-06-30_06-30-27');
aveTreat_GeneXinsteadOfSS = strcat(masterFolder, 'inf0p0023_lat2p16_base_Post2011_GeneXSS\2012-06-30_10-13-33');


%%%% GRAB INPUT MATRACIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%base case stuff
cd(aveTreat_fullCat4_empUptake);
aveTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
aveTreat_fullCat4_empUptake_incidence = dlmread('incidenceMatrix_forMakingWHOComparison.csv', ',' ,1,0);

aveTreat_fullCat4_empUptake_ageArray{1} = dlmread('HealthOutcomes_age0_20Only.csv', ',' ,1,0);
aveTreat_fullCat4_empUptake_ageArray{2} = dlmread('HealthOutcomes_age21_60Only.csv', ',' ,1,0);
aveTreat_fullCat4_empUptake_ageArray{3} = dlmread('HealthOutcomes_age61_100Only.csv', ',' ,1,0);

aveTreat_fullCat4_empUptake_smokeArray{1} = dlmread('HealthOutcomes_nonSmokingOnly.csv', ',' ,1,0);
aveTreat_fullCat4_empUptake_smokeArray{2} = dlmread('HealthOutcomes_smokingOnly.csv', ',' ,1,0);

%no treat
cd(noTreat);
noTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

noTreat_fullCat4_empUptake_ageArray{1}  = dlmread('HealthOutcomes_age0_20Only.csv', ',' ,1,0);   %age 0_20
noTreat_fullCat4_empUptake_ageArray{2}  = dlmread('HealthOutcomes_age21_60Only.csv', ',' ,1,0);   %age 21_60
noTreat_fullCat4_empUptake_ageArray{3}  = dlmread('HealthOutcomes_age61_100Only.csv', ',' ,1,0);  %age 61_100

noTreat_fullCat4_empUptake_smokeArray{1} = dlmread('HealthOutcomes_nonSmokingOnly.csv', ',' ,1,0);     %nonSmoke
noTreat_fullCat4_empUptake_smokeArray{2} = dlmread('HealthOutcomes_smokingOnly.csv', ',' ,1,0);    %smoking

% %uptake stuff
% cd(aveTreat_noCat4_ultraUptake)
% aveTreat_noCat4_ultraUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
%
cd(aveTreat_noCat4_empUptake)
aveTreat_noCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
%
% %MDR evolved vs. transmitted
cd(aveTreat_fullCat4_empUptake);
baseCase_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);

cd(noEvolve)
noEvolve_2011start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
noEvolve_MDRevolvedTrans = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);

%
% %MDR and treatment quality
cd(perfectTreat_fullCat4_empUptake);
perfectTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

cd(bestTreat_fullCat4_empUptake);
bestTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

cd(worstTreat_fullCat4_empUptake);
worstTreat_fullCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

%
% %MDR under different DOTS+ scenarios
cd(bestTreat_perfectCat4_empUptake);
bestTreat_perfectCat4_empUptake_outcomes  = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(aveTreat_perfectCat4_empUptake);
aveTreat_perfectCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(bestTreat_noCat4_empUptake);
cbestTreat_noCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(perfectTreat_perfectCat4_empUptake)
perfectTreat_perfectCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(perfectTreat_noCat4_empUptake)
perfectTreat_noCat4_empUptake_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

%%%%half dots plus stuff
cd(halfDotsPlus);
halfDotsPlus_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(halfDotsPlusNoEVOLVE_2013start)
halfDotsPlus_noEvolve2013_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(halfDotsPlusNoEVOLVE_2019start)
halfDotsPlus_noEvolve2019_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(halfDotsPlusNoEVOLVE_2027start)
halfDotsPlus_noEvolve2027_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

% DST scenarios
cd(aveTreat_allGeneXscenarios);
aveTreat_allGeneXscenarios_outcomes  = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(aveTreat_GeneXinitialDST);
aveTreat_GeneXinitialDST_outcomes  = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(aveTreat_GeneXinsteadOfLJ);
aveTreat_GeneXinsteadOfLJ_outcomes  = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(aveTreat_GeneXinsteadOfSS);
aveTreat_GeneXinsteadOfSS_outcomes  = dlmread('HealthOutcomes.csv', ',' ,1,0);


% %last noEvolve
cd(aveTreat_noEvolve_2019start)
aveTreat_noEvolve_2019start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(aveTreat_noEvolve_2027start)
aveTreat_noEvolve_2027start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(bestTreat_2019start)
bestTreat_2019start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);
cd(bestTreat_2027start)
bestTreat_2027start_outcomes = dlmread('HealthOutcomes.csv', ',' ,1,0);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTTING
if startYear == 26
    startCalYear = 1996;
    startYroffset = 0; %for the MDR stuff -- the numMDRevolvedOrTransmitted file starts at 1996 and is monthly
elseif startYear == 29
    startCalYear = 2011;
    startYroffset = 180; %for the MDR stuff -- the numMDRevolvedOrTransmitted file starts at 1996 and is monthly so if you want to start at 2011 gotta start at row 180
end
years = [startCalYear:5:2036]';

%%proprotions by age/smoking graph

for typeTB = 1:4  %loop over latDS, actDS, latMDR, actMDR
    
    if typeTB == 1
        healthOutcomesCol = 7; %do latent DS
        TBtypeTitleStr = 'Latent DS';
        ylimMat = [0 0.9];
    elseif typeTB == 2
        healthOutcomesCol = 9; %active DS
        TBtypeTitleStr = 'Active DS';
        ylimMat = [0 0.012];
    elseif typeTB == 3
        healthOutcomesCol = 8; %latent MDR
        TBtypeTitleStr = 'Latent MDR';
        ylimMat = [0 0.025];
    elseif typeTB == 4
        healthOutcomesCol = 10; %act MDR
        TBtypeTitleStr = 'Active MDR';
        ylimMat = [0 0.0008];
    end
    
    for demogType = 1:3  %loop over demographics (smoking, age etc)
        %set demographic type specifics
        if demogType == 1
            demogTypeTitleStr = 'Smoking Status';
            numTypes = 2; %either nonsmoking or smoking
            for i = 1:numTypes
                noTreatDemogGrp{i} = noTreat_fullCat4_empUptake_smokeArray{i};
                aveTreatDemogGrp{i} = aveTreat_fullCat4_empUptake_smokeArray{i};
            end
        elseif demogType == 2
            demogTypeTitleStr = 'Age Group';
            numTypes = 3; %age under20, 21-60, or 61+
            for i = 1:numTypes
                noTreatDemogGrp{i} = noTreat_fullCat4_empUptake_ageArray{i};
                aveTreatDemogGrp{i} = aveTreat_fullCat4_empUptake_ageArray{i};
            end
        elseif demogType == 3
            demogTypeTitleStr = 'Overall';
            numTypes = 1; %overall DS TB graph
            for i = 1:numTypes
                noTreatDemogGrp{i} = noTreat_fullCat4_empUptake_outcomes;
                aveTreatDemogGrp{i} = aveTreat_fullCat4_empUptake_outcomes;
            end
            
        end
        
        %titles and labels
        TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
        demogTypeStr = sscanf(demogTypeTitleStr,'%s'); %demog with no spaces
        titleStr = strcat('Proportion of Population with', TBtypeTitleStr,' Over Time, By ', demogTypeTitleStr);
        printStr = strcat('Prev','_', TBtype, '_', demogTypeStr);
        
        %color and line types
        lineStyle_NoTreat = '--';
        lineStyle_Treat = '-';
        
        if demogType == 1
            treatColor{1} = [0.65 0.50 0.25];
            treatColor{2} = [0.65 0.50 1];
        elseif demogType == 2
            treatColor{1} = 'r';
            treatColor{2} = 'm';
            treatColor{3} = 'b';
        elseif demogType == 3
            treatColor{1} = [0 0 0];
        end
        
        %get data for plotting
        for y = 1:numTypes
            propAveTreat{y} = aveTreatDemogGrp{y}(startYear:34, healthOutcomesCol);  %latent DS
            propNoTreat{y} = noTreatDemogGrp{y}(startYear:34, healthOutcomesCol);  %latent DS
        end
        
        %plot everything and print
        for y = 1:numTypes
            cd(pptGraphsFold);
            plot(years, propNoTreat{y}, 'Color', treatColor{y}, 'Linestyle', lineStyle_NoTreat);
            hold on
        end
        ylim(ylimMat);
        cd(plotFixerFold);
        plotfixer;
        cd(pptGraphsFold);
        plotName = strcat('noTreat', printStr);
        print('-dpng',plotResolution,[pptGraphsFold '/' plotName]);
        close all
        
        %generate the same graph again except with the aveTreat overlaid
        for y = 1:numTypes
            plot(years, propNoTreat{y},'Color', treatColor{y}, 'Linestyle', lineStyle_NoTreat);
            hold on
            plot(years, propAveTreat{y},'Color', treatColor{y}, 'Linestyle', lineStyle_Treat);
            hold on
        end
        ylim(ylimMat);
        cd(plotFixerFold);
        plotfixer;
        plotName = strcat('no_andAveTreat', printStr);
        print('-dpng',plotResolution,[pptGraphsFold '/' plotName]);
        close all
        
    end %end demog type
end  %end TB type


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%transmission vs evolved graph

months = [1996:1/12:(2036+11/12)]';

fileName{1} = baseCase_MDRevolvedTrans; caseStr{1} = 'base';
fileName{2} = noEvolve_MDRevolvedTrans; caseStr{2} = 'noEvolve';

extraOverlayGraph = 0;
yearly_MDRevolvedTrans = zeros(floor(((2036-startCalYear+1)*12)/60),2);  %because numNewMDR_evolvedOrTrans_postBurnIn is in months from 1996 to 2046, so need to stop at 492, the end of 2036
%sum to get five year bins
for caseNum = 1:2
    
    for i = 0:(floor(((2036-startCalYear+1)*12)/60)-1)
        yearly_MDRevolvedTrans(i+1,:) = sum( fileName{caseNum}( (60*i+1+startYroffset):60*i+60+startYroffset, 1:2) );
    end
    
    totNumPpl = sum(aveTreat_fullCat4_empUptake_outcomes(startYear:33,1:6)')';
    totNewMDRcases = sum(yearly_MDRevolvedTrans')';
    percEvolved = yearly_MDRevolvedTrans(:,1)./totNewMDRcases(:,1);
    percTransmitted = yearly_MDRevolvedTrans(:,2)./totNewMDRcases(:,1);
    
    toPlot = [yearly_MDRevolvedTrans(:,1)./totNumPpl, yearly_MDRevolvedTrans(:,2)./totNumPpl];
    
    %     plot(years, yearly_MDRevolvedTrans(:,1)./totNumPpl, 'r');
    %     %plot(years, percEvolved, 'r');
    %     hold on;
    %     plot(years, totNewMDRcases./totNumPpl, 'b');
    %
    
    bar(years(1:size(years,1)-1,:), toPlot(1:size(toPlot,1), :), 'stack');
    xlim([2008 2035]);
    ylim([0 0.00175]);
    
    evolvedTransmittedGraphs{caseNum} = toPlot(1:size(toPlot,1), :);
    P=findobj(gca,'type','patch');
    if caseNum == 1
        C=['m','b','r','g']; % make a colors list
        for n=1:length(P)
            set(P(n),'facecolor',C(n));
        end
    elseif caseNum ==2
        C=['m','g','r','b']; % make a colors list
        for n=1:length(P)
            set(P(n),'facecolor',C(n));
        end
    end
    
    %plot(years, percTransmitted, 'b');
    % xlabel('Year');
    % ylabel('Proportion of New MDR Cases');
    % title('Proportion of New MDR Cases Evolved vs. Transmitted');
    % legend('Evolved', 'Transmitted');
    cd(plotFixerFold);
    plotfixer;
    fileStr = strcat('MDRevolvedTransmitted_', caseStr{caseNum});
    print('-dpng',plotResolution,[pptGraphsFold '/' fileStr]);
    close all
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%MDR and treatment quality
for typeTB = 1:2  %loop over latMDR, actMDR
    
    if typeTB == 1
        healthOutcomesCol = 8; %latent MDR
        TBtypeTitleStr = 'Latent MDR';
        ylimVec = [0 0.03];
    elseif typeTB == 2
        healthOutcomesCol = 10; %act MDR
        TBtypeTitleStr = 'Active MDR';
        ylimVec = [0 0.0008];
    end
    TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
    titleStr = strcat('Proportion of Population with ', TBtype, ' , by DOTS Quality');
    
    aveTreat = aveTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    worstTreat = worstTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    bestTreat = bestTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    perfectTreat = perfectTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    
    for graphLayer = 1:4
        plot(years, aveTreat, 'Color', [0.1 0.5 0.9]);
        layerNum = '1';
        if graphLayer > 1
            hold on;
            plot(years, worstTreat, 'b');
            layerNum = '2';
        end
        if graphLayer > 2
            hold on;
            plot(years, bestTreat, 'r');
            layerNum = '3';
        end
        if graphLayer > 3
            hold on;
            plot(years, perfectTreat, 'Color',[0.9 0.5 0.1]);
            layerNum = '4';
        end
        ylim(ylimVec);
        cd(plotFixerFold);
        plotfixer;
        printStr = strcat('treatmentQuality_', TBtype, '_', layerNum);
        print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
        close all
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%MDR under different DOTS+ scenarios


for graphRun = 1:3  %make three graphs so can do appear on ppt
    for typeTB = 1:2  %loop over latMDR, actMDR
        
        if typeTB == 1
            healthOutcomesCol = 8; %latent MDR
            TBtypeTitleStr = 'Latent MDR';
        elseif typeTB == 2
            healthOutcomesCol = 10; %act MDR
            TBtypeTitleStr = 'Active MDR';
        end
        TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
        
        perfectPerf = perfectTreat_perfectCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        perfectNormalCat4 = perfectTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        perfectNo = perfectTreat_noCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        
        avePerf = aveTreat_perfectCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        aveNormalCat4 = aveTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        aveNo = aveTreat_noCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        
        if  graphRun == 3
            printStr = strcat('MDRcontrol_', TBtype, '_wPerfCat4');
        elseif graphRun == 2
            printStr = strcat('MDRcontrol_', TBtype, '_wNormalCat4');
        elseif graphRun == 1
            printStr = strcat('MDRcontrol_', TBtype, '_noCatIV');
        end
        
        if graphRun > 2
            plot(years, perfectPerf, 'Color', [0.9 0.5 0.1],'Linestyle', '-.');
            hold on;
            plot(years, avePerf, 'Color', [0.1 0.5 0.9], 'Linestyle', '-.');
            hold on;
        end
        if graphRun > 1
            plot(years, perfectNormalCat4, 'Color', [0.9 0.5 0.1],'Linestyle', '--');
            hold on;
            plot(years, aveNormalCat4, 'Color', [0.1 0.5 0.9], 'Linestyle', '--');
            hold on;
        end
        plot(years, perfectNo, 'Color', [0.9 0.5 0.1]);
        hold on;
        plot(years, aveNo, 'Color', [0.1 0.5 0.9]);
        
        if typeTB == 1
            ylim([0 0.04]);
        elseif typeTB == 2
            ylim([0 0.0005]);
        end
        
        cd(plotFixerFold);
        plotfixer;
        print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
        close all
        
    end %end TBtype
    
end %end graphRun


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DST graphs
for graphRun = 1:5  %make 5 graphs so can do appear on ppt
    for typeTB = 1:4  %loop over latDS, actDS, latMDR, actMDR
        
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
        TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
        
        if  graphRun == 5
            printStr = strcat('DST_', TBtype, '_aveTreat_allGeneXscenarios'); %'_aveTreat_allGeneXscenarios'
        elseif graphRun == 4
            printStr = strcat('DST_', TBtype, '_aveTreat_initialDST');  %'_aveTreat_GeneXinitialDST'
        elseif graphRun == 3
            printStr = strcat('DST_', TBtype, '_aveTreat_SSisPerfect');  %'_aveTreat_GeneXinsteadOfSS'
        elseif graphRun == 2
            printStr = strcat('DST_', TBtype, '_aveTreat_fastDST');  %'_aveTreat_GeneXinsteadOfLJ'
        elseif graphRun == 1
            printStr = strcat('DST_', TBtype, '_baseOnly');
        end
        
        allGeneX_run = aveTreat_allGeneXscenarios_outcomes(startYear:34,healthOutcomesCol);
        aveTreat_GeneXinitialDST_run = aveTreat_GeneXinitialDST_outcomes(startYear:34,healthOutcomesCol);
        aveTreat_GeneXinsteadOfSS_run = aveTreat_GeneXinsteadOfSS_outcomes(startYear:34,healthOutcomesCol);
        aveTreat_GeneXinsteadOfLJ_run = aveTreat_GeneXinsteadOfLJ_outcomes(startYear:34,healthOutcomesCol);
        basecaseRun = aveTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
        
        if graphRun > 4
            plot(years, allGeneX_run, 'Color', [0.7 0 1]);
            hold on;
        end
        if graphRun > 3
            plot(years, aveTreat_GeneXinitialDST_run, 'Color', [0.7 0.3 0.8]);
            hold on;
        end
        if graphRun > 2
            plot(years, aveTreat_GeneXinsteadOfSS_run, 'Color', [0.7 0.6 0.6]);
            hold on;
        end
        if graphRun > 1
            plot(years, aveTreat_GeneXinsteadOfLJ_run, 'Color', [0.7 0.9 0.4]);
            hold on;
        end
        plot(years, basecaseRun, 'Color', [0.1 0.6 0.2]);
        
        cd(plotFixerFold);
        plotfixer;
        print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
        close all
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WHO comparison graphs, reformatted

timeAxis = [1986:5:2011];
propActDS = aveTreat_fullCat4_empUptake_outcomes(24:29,9); %starting in 1986, actDS

%do prevalence by rate (ie, total active out of 100,000) and put
%the WHO numbers on too for comparison
dataYear = [1990	; 1991	; 1992	; 1993	; 1994	;  1995	; 1996	; 1997	; 1998	; 1999	; 2000	; 2001	; 2002	; 2003	; 2004	; 2005	; 2006	; 2007	; 2008	; 2009	; 2010	];
actTBPrev = [459 ; 460 ; 460	; 461	; 462	; 462	; 463	; 464	; 464	; 465	; 466	; 466	; 436	; 409	; 383	; 358	; 335	; 314	; 294	; 275	; 256];
distToUpperLimit = [56	; 56	; 56	; 56	; 56	; 57	; 56	; 56	; 57	; 56	; 56	; 57	; 62	; 66	; 71	; 78	; 84	; 91	; 99	; 107	; 117];
distToLowerLimit = [52	; 53	; 52	; 53	; 53	; 53	; 53	; 53	; 53	; 53	; 54	; 53	; 57	; 62	; 66	; 70	; 74	; 80	; 85	; 90	; 95];

%make simulation prevalences into rates out of 100000
for graphLayer = 1:2
    titleStrSuffix = 'WHOonly';
    %errorbar(dataYear, actTBPrev,distToLowerLimit,distToUpperLimit);
    plot(dataYear, actTBPrev + distToUpperLimit, 'b-.')
    hold on;
    plot(dataYear, actTBPrev - distToLowerLimit, 'b-.')
    hold on;
    plot(dataYear, actTBPrev, 'b')
    
    if graphLayer == 2
        hold on;
        toPlot = (propActDS) *100000;
        p = plot(timeAxis,toPlot);
        set(p,'Color','red','LineWidth',1);
        titleStrSuffix = 'withPrevLine';
    end
    %     title('Total active TB, out of 100000, with WHO estimates');
    %     xlabel('Year');
    %     ylabel('Prevalence Rate');
    ylim([0 600])
    xlim([1980 2020]);
    cd(plotFixerFold);
    plotfixer;
    titleString = strcat('TBprevWHOcomparison_', titleStrSuffix, '.png');
    print('-dpng',plotResolution,[pptGraphsFold '/' titleString]);
    close all
end


%
timeAxis = [1986:1:2011];
incidence_monthly = aveTreat_fullCat4_empUptake_incidence(1440:1751, :);  %first column is incidence, in months, and column 3 is total pop.  1440 to 1751 is 1986 to 2011

%sum into yearly bins
incidence_yearly = zeros( (2011-1986+1) , size(incidence_monthly,2) );
for binNum = 0:(2011-1986)
    incidence_yearly(binNum+1,:) = sum( incidence_monthly(binNum*12 + 1 : binNum*12 + 12 ,:) );
    incidence_yearly(binNum+1,3) = incidence_monthly(binNum*12 + 7,3) ;
end

toPlot = (  incidence_yearly(:,1) ./ incidence_yearly(:,3) )* 100000;

%do incidence by rate (ie, total active out of 100,000) and put
dataTime = [1990	; 1991	; 1992	; 1993	; 1994	;  1995	; 1996	; 1997	; 1998	; 1999	; 2000	; 2001	; 2002	; 2003	; 2004	; 2005	; 2006	; 2007	; 2008	; 2009	; 2010	];
actTBinci = [216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 216	; 215	; 214	; 212	; 209	; 205	; 201	; 196	; 190	; 185 ];
distToUpperLimit = [39	; 37	; 35	; 33	; 32	; 30	; 29	; 27	; 26	; 25	; 24	; 23	; 23	; 23	; 22	; 22	; 22	; 21	; 21	; 21	; 20	];
distToLowerLimit = [35	;  33	; 32	; 30	; 29	; 27	; 26	; 25	; 24	; 23	; 22	; 22	; 22	; 22	; 22	;21	    ; 21	; 20	; 20	; 19	; 18	];

%make simulation incidences into rates out of 100000 and plot
for graphLayer = 1:2
    titleStrSuffix = 'WHOonly';
    %errorbar(dataTime, actTBinci,distToLowerLimit,distToUpperLimit);
    plot(dataTime,actTBinci+distToUpperLimit, 'b-.' );
    hold on;
    plot(dataTime,actTBinci-distToLowerLimit, 'b-.' );
    hold on;
    plot(dataTime,actTBinci, 'b' );
    if graphLayer == 2
        hold on;
        p = plot(timeAxis,toPlot);
        titleStrSuffix = 'withIncLine';
        set(p,'Color','red','LineWidth',1);
    end
    ylim([ 0 260 ]);
    xlim([1980 2020]);
    cd(plotFixerFold);
    plotfixer;
    titleString = strcat('TBincidence_WHOcomparison', titleStrSuffix, '.png');
    print('-dpng',plotResolution,[pptGraphsFold '/' titleString]);
    close all
end



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
    
    %   if graphRun > 1
    %       plot(years, fastDST, 'Color', [0.7 0.5 0.4]);
    %       hold on;
    %   end
    for graphLayer = 1:4
        layerStr = '1';
        if graphLayer == 2
            layerStr = '2';
        elseif graphLayer == 3
            layerStr = '3';
        elseif graphLayer == 4
            layerStr = '4';
        end
        
        
        if graphLayer >= 4
            plot(years,  start2027 , 'Color', [0.5 0 0.5]);
            hold on;
        end
        if graphLayer >= 3
            plot(years, start2019, 'Color', [0.9 0.4 0.2]);
            hold on;
        end
        if graphLayer >= 2
            plot(years, start2011 , 'Color',  [0.9 0.5 0.1]);
            hold on;
        end
        plot(years, base, 'Color',  [0.1 0.5 0.9]);
        
        %
        %             if typeTB == 1
        %                 ylim([0 0.02]);
        %             elseif typeTB == 2
        %                 ylim([0 0.003]);
        %             end
        ylim([0 0.0003])
        cd(plotFixerFold);
        plotfixer;
        printStr = strcat('delayInNoEvolve_', TBtypeTitleStr, '_', layerStr);
        print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
        close all
    end
    
    
end



%The no evolve graphs, except using bestTreat instead
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
    start2011 =  bestTreat_fullCat4_empUptake_outcomes(startYear:34, healthOutcomesCol);
    start2019 = bestTreat_2019start_outcomes(startYear:34, healthOutcomesCol);
    start2027 = bestTreat_2027start_outcomes(startYear:34, healthOutcomesCol);
    
    %   if graphRun > 1
    %       plot(years, fastDST, 'Color', [0.7 0.5 0.4]);
    %       hold on;
    %   end
    for graphLayer = 1:4
        layerStr = '1';
        if graphLayer == 2
            layerStr = '2';
        elseif graphLayer == 3
            layerStr = '3';
        elseif graphLayer == 4
            layerStr = '4';
        end
        
        
        if graphLayer >= 4
            plot(years,  start2027 , 'Color', [0.5 0 0.5]);
            hold on;
        end
        if graphLayer >= 3
            plot(years, start2019, 'Color', [0.9 0.4 0.2]);
            hold on;
        end
        if graphLayer >= 2
            plot(years, start2011 , 'Color',  [0.9 0.5 0.1]);
            hold on;
        end
        plot(years, base, 'Color',  [0.1 0.5 0.9]);
        
        %
        %             if typeTB == 1
        %                 ylim([0 0.02]);
        %             elseif typeTB == 2
        %                 ylim([0 0.003]);
        %             end
        ylim([0 0.0003])
        cd(plotFixerFold);
        plotfixer;
        printStr = strcat('delayInBestTreat_', TBtypeTitleStr, '_', layerStr);
        print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
        close all
    end
    
    
end




%%%MDR under different DOTS+, bestTreat, and noEvolve scenarios

%cases:
%base
%base without DOTS+
%best2011
%best2019
%best2027
%noEvolve2011
%noEvolve2019
%noEvolve2027

%halfDotsPlus
%halfDotsPlus_noEvolve2013
%halfDotsPlus_noEvolve2019
%halfDotsPlus_noEvolve2027


for typeTB = 1:2  %loop over latMDR, actMDR
    
    if typeTB == 1
        healthOutcomesCol = 8; %latent MDR
        TBtypeTitleStr = 'Latent MDR';
    elseif typeTB == 2
        healthOutcomesCol = 10; %act MDR
        TBtypeTitleStr = 'Active MDR';
    end
    TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
    
    
    aveNormalCat4 = aveTreat_fullCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    aveNoCat4 = aveTreat_noCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    avePerfectCat4 = aveTreat_perfectCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
    
    bestStart2011 =  bestTreat_fullCat4_empUptake_outcomes(startYear:34, healthOutcomesCol);
    bestStart2019 = bestTreat_2019start_outcomes(startYear:34, healthOutcomesCol);
    bestStart2027 = bestTreat_2027start_outcomes(startYear:34, healthOutcomesCol);
    
    noEvolveStart2011 =  noEvolve_2011start_outcomes(startYear:34, healthOutcomesCol);
    noEvolveStart2019 = aveTreat_noEvolve_2019start_outcomes(startYear:34, healthOutcomesCol);
    noEvolveStart2027 = aveTreat_noEvolve_2027start_outcomes(startYear:34, healthOutcomesCol);
    
    halfDotsPlus  = halfDotsPlus_outcomes(startYear:34, healthOutcomesCol);
    halfDotsPlus_noEvolve2013 = halfDotsPlus_noEvolve2013_outcomes(startYear:34, healthOutcomesCol);
    halfDotsPlus_noEvolve2019 = halfDotsPlus_noEvolve2019_outcomes(startYear:34, healthOutcomesCol); 
    halfDotsPlus_noEvolve2027 = halfDotsPlus_noEvolve2027_outcomes(startYear:34, healthOutcomesCol); 

    
    
    printStr = strcat('MDRsummary_', TBtype);
    
    plot(years, aveNoCat4, 'Color', [0.5 0 0.5]);
    hold on;
    plot(years, halfDotsPlus, 'Color', [0.3 0.3 0.6]);
    hold on;
    plot(years, aveNormalCat4, 'Color', [0.1 0.5 0.9]);
    hold on;
    plot(years, avePerfectCat4, 'Color', [0.1 0.6 0.2]);
    
    plot(years, bestStart2027, 'Color', 'r','Linestyle', '-.');
    hold on;
    plot(years, bestStart2019, 'Color', 'r', 'Linestyle', '--');
    hold on;
    plot(years, bestStart2011, 'Color', 'r');
    hold on;
    
    plot(years, halfDotsPlus_noEvolve2027, 'Color', [0.9 0.9 0.3],'Linestyle', '-.');
    hold on;
    plot(years, halfDotsPlus_noEvolve2019, 'Color', [0.9 0.9 0.3], 'Linestyle', '--');
    hold on;
    plot(years, halfDotsPlus_noEvolve2013, 'Color', [0.9 0.9 0.3]);
    hold on;
    
    plot(years, noEvolveStart2027, 'Color', [0.9 0.5 0.1],'Linestyle', '-.');
    hold on;
    plot(years, noEvolveStart2019, 'Color', [0.9 0.5 0.1], 'Linestyle', '--');
    hold on;
    plot(years, noEvolveStart2011, 'Color', [0.9 0.5 0.1]);
    hold on;
    
    h_leg = legend('No DOTS+', 'Half DOTS+', 'DOTS+', 'Ideal DOTS+', 'Best Quality 2027','Best Quality 2019','Best Quality 2012','Half DOTS+, No Treatment Gen. MDR 2027','Half DOTS+, No Treatment Gen MDR 2019','Half DOTS+, No Treatment Gen MDR 2012', 'No Treatment Generated MDR 2027','No Treatment Generated MDR 2019','No Treatment Generated MDR 2012', 'Location', 'NorthWest' );
    
    %       if typeTB == 1
    %            ylim([0 0.04]);
    %        elseif typeTB == 2
    %            ylim([0 0.0005]);
    %        end
    
    cd(plotFixerFold);
    plotfixer;
    set(h_leg,'FontSize',11);
    print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
    close all
    
end %end TBtype










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBUGGING AND OTHER SCRAP CODE
%
% THIS WOULD GO AFTER THE EVOLVED VS TRANSMITTED SECTION
%
% interMedGraph = evolvedTransmittedGraphs{1} - evolvedTransmittedGraphs{2};
% interMedGraph = [interMedGraph, evolvedTransmittedGraphs{2}];
% bar(years(1:size(years,1)-1,:), interMedGraph, 'stack');
% P=findobj(gca,'type','patch');
% C=['m','g','r','c']; % make a colors list
% for n=1:length(P)
%     set(P(n),'facecolor',C(n));
% end
% xlim([2008 2035]);
% ylim([0 0.00175]);
% cd(plotFixerFold);
% plotfixer;
% fileStr = strcat('MDRevolvedTransmitted_', 'diff_if_noEvolve');
% print('-dpng',plotResolution,[pptGraphsFold '/' fileStr]);
% close all
%
%




%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Higher uptake means more MDR
% for typeTB = 1:4  %loop over latDS, actDS, latMDR, actMDR
%
%
% if typeTB == 1
%         healthOutcomesCol = 7; %do latent DS
%         TBtypeTitleStr = 'Latent DS';
%     elseif typeTB == 2
%         healthOutcomesCol = 9; %active DS
%         TBtypeTitleStr = 'Active DS';
%     elseif typeTB == 3
%         healthOutcomesCol = 8; %latent MDR
%         TBtypeTitleStr = 'Latent MDR';
%     elseif typeTB == 4
%         healthOutcomesCol = 10; %act MDR
%         TBtypeTitleStr = 'Active MDR';
%     end
%
%
%     TBtype = sscanf(TBtypeTitleStr,'%s');%TB type with no spaces
%     titleStr = strcat(TBtype, ' Prevalence, by Uptake Level');
%     printStr = strcat('uptakeComp_', TBtype);
%
%     ultraUptake = aveTreat_noCat4_ultraUptake_outcomes(startYear:34,healthOutcomesCol);
%     %womUptake = aveTreat_noCat4_womUptake_outcomes(startYear:34,healthOutcomesCol);
%     empUptake = aveTreat_noCat4_empUptake_outcomes(startYear:34,healthOutcomesCol);
%
%     plot(years, ultraUptake, 'Color', [0.60 0.50 0.90]);
%     hold on;
%     %         plot(years, womUptake, 'm');
%     %         hold on;
%     plot(years, empUptake, 'Color', [0.55 0.75 0]);
%
%     cd(plotFixerFold);
%     plotfixer;
%     print('-dpng',plotResolution,[pptGraphsFold '/' printStr]);
%     close all
% end



