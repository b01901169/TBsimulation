function smokingPrevPlotter(stateMat, numCurrentPpl, timePeriod, Page, Psex, Psmoking, plotResolution, folderName, subset)
%function call is
% smokingPrevPlotter(stateMat,numCurrentPpl, subset)
% if subset is 0, no subset
% if subset is 1, subset is sex (male and female graphs generated)

if subset == 0
    %no subset
    tots = zeros(1,21);
    smokersTots = zeros(1,21);
    for indivs = 1:numCurrentPpl
        ageGrp = floor(double(stateMat(indivs,Page)) / 5) + 1;
        tots(1,ageGrp) = tots(1,ageGrp) + 1;
        if stateMat(indivs,Psmoking) == 1
            smokersTots(1,ageGrp) = smokersTots(1,ageGrp) + 1;
        end
    end

    subplot(2,1,1), plot([0:5:100], tots, 'r');
    hold on;
    plot([0:5:100], smokersTots,'b');
    axis([0 100 0 350000])
    legend('total number in age grp','smokers in age grp');
    title('Smoking Frequency Over Age');

    toPlot = smokersTots./tots;
    subplot(2,1,2), plot([0:5:100],toPlot);
    axis([0 100 0 0.25])
    title('Smoking Proportion Over Age');
    smokingplottitle = sprintf('smoking%03d',timePeriod);
    print('-dpng',plotResolution,[folderName '/' smokingplottitle '.png']);
    close all
    
    tablePrinter('numSmokers by 5yr age groups', smokersTots, 'smokersTots', folderName)
    
elseif subset == 1
    %subsetting on sex
    totsMale = zeros(1,21);
    smokersTotsMale = zeros(1,21);
    totsFemale = zeros(1,21);
    smokersTotsFemale = zeros(1,21);
	totPredMaleSmokers = zeros(1,21);
	totPredFemaleSmokers = zeros(1,21);

    for indivs = 1:numCurrentPpl
        if stateMat(indivs,Psex) == 1
            ageGrp = floor(double(stateMat(indivs,Page)) / 5) + 1;
            totsMale(1,ageGrp) = totsMale(1,ageGrp) + 1;
            if stateMat(indivs,Psmoking) == 1
                smokersTotsMale(1,ageGrp) = smokersTotsMale(1,ageGrp) + 1;
            end
        elseif stateMat(indivs,Psex) == 2
            ageGrp = floor(double(stateMat(indivs,Page)) / 5) + 1;
            totsFemale(1,ageGrp) = totsFemale(1,ageGrp) + 1;
            if stateMat(indivs,Psmoking) == 1
                smokersTotsFemale(1,ageGrp) = smokersTotsFemale(1,ageGrp) + 1;
            end
        end
    end
    
    %generated predicted num of smokers
%     
%     smokingPrev =[...
%         0	0
%         0	0
%         0	0
%         0.06	0.02
%         0.12	0.025
%         0.22	0.03
%         0.28	0.035
%         0.35	0.03
%         0.36	0.035
%         0.37	0.04
%         0.4     0.04
%         0.38	0.05
%         0.38	0.06
%         0.34	0.055
%         0	0
%         0	0
%         0	0
%         0	0
%         0	0
%         0	0
%         0	0
%         ];  %from Jha NEJM (read off figure 1) for data in 2002ish.  we want to match this, generally.
    

smokingPrev = [... 
0.0000000	0.0000000
0.0000000	0.0000000
0.0000000	0.0000000
0.0540131	0.0025893
0.1617864	0.0068130
0.2213199	0.0141294
0.2808533	0.0214457
0.3190453	0.0301124
0.3572372	0.0387790
0.3712019	0.0489459
0.3851665	0.0591128
0.3597671	0.0727734
0.3343677	0.0864339
0.3258945	0.0958011
0.3174212	0.1051682
0.2916656	0.0959240
0.2574367	0.0960470
0.2627954	0.1343809
0.2681540	0.1727147
0.1457413	0.1142369
0.0233286	0.0557590
];  %replaced the with the GATS numbers that Jeremy sent me.  95% CI are made from the binomial exact confidence interval eqn.

smokingPrevLowerBound = [...
0.0000000	0.0000000
0.0000000	0.0000000
0.0000000	0.0000000
0.0456463	0.0007698
0.1533054	0.0051675
0.2124520	0.0118647
0.2715986	0.0185619
0.3088206	0.0262658
0.3460427	0.0339697
0.3578168	0.0425799
0.3695909	0.0511901
0.3425183	0.0631094
0.3154456	0.0750287
0.3028634	0.0803446
0.2902811	0.0856604
0.2547692	0.0712972
0.2066751	0.0622499
0.1783358	0.0688411
0.1499965	0.0754323
0.0398842	-0.0022042
-0.0702280	-0.0798407
];

smokingPrevUpperBound = [...
0.0000000	0.0000000
0.0000000	0.0000000
0.0000000	0.0000000
0.0623799	0.0044088
0.1702674	0.0084585
0.2301877	0.0163940
0.2901080	0.0243295
0.3292699	0.0339589
0.3684317	0.0435883
0.3845869	0.0553119
0.4007421	0.0670355
0.3770159	0.0824373
0.3532898	0.0978391
0.3489255	0.1112575
0.3445613	0.1246760
0.3285619	0.1205508
0.3081983	0.1298441
0.3472549	0.1999206
0.3863115	0.2699971
0.2515984	0.2306779
0.1168852	0.1913587
];


    totPredMaleSmokers = totsMale .* smokingPrev(:,1)';
    totPredFemaleSmokers = totsFemale .* smokingPrev(:,2)';

    %generate titles by sex
    for subsetValue = 1:2
        if subsetValue == 1
            title1 = 'Male Smoking Frequency Over Age';
            title2 = 'Male Smoking Proportion Over Age';
            smokingplottitle = sprintf('Male_smoking%03d',timePeriod);
            tots = totsMale;
            smokersTots = smokersTotsMale;
            predSmokersTots = totPredMaleSmokers;
            toPlotPred = smokingPrev(:,1);
            toPlotPredUpperBound = smokingPrevUpperBound(:,1);
            toPlotPredLowerBound = smokingPrevLowerBound(:,1);
        elseif subsetValue == 2
            title1 = 'Female Smoking Frequency Over Age';
            title2 = 'Female Smoking Proportion Over Age';
            smokingplottitle = sprintf('Female_smoking%03d',timePeriod);
            tots = totsFemale;
            smokersTots = smokersTotsFemale;
            predSmokersTots = totPredFemaleSmokers;
            toPlotPred = smokingPrev(:,2);
            toPlotPredUpperBound = smokingPrevUpperBound(:,2);
            toPlotPredLowerBound = smokingPrevLowerBound(:,2);
        end

        %plot the graphs
        subplot(2,1,1), plot([0:5:100], tots, 'r');
        hold on;
        plot([0:5:100], smokersTots,'b');
        hold on;
        plot([0:5:100],predSmokersTots,'g');
        axis([0 100 0 500000])
        legend('total number in age grp','smokers in age grp', 'smokers if not use churn');
        title(title1);

        toPlot = smokersTots./tots;
        subplot(2,1,2), plot([0:5:100],toPlot,'b');
        hold on;
        plot([0:5:100], toPlotPred,'g');
        hold on;
        plot([0:5:100], toPlotPredUpperBound,'g-.');
        hold on;
        plot([0:5:100], toPlotPredLowerBound,'g-.');
        axis([0 100 0 1])
        title(title2);
        print('-dpng',plotResolution,[folderName '/' smokingplottitle '.png']);
        close all
        
    	tablePrinter('num Smokers by 5yr age groups', smokersTots, smokingplottitle, folderName)

    end
end
