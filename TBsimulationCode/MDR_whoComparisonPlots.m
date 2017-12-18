function MDRoutVals = MDR_whoComparisonPlots(folderName)
makeGraphs = 0;  %graph are all hard coded


% 
% folderName{1,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_1997\b01_2013-02-13_18-12-58';
% folderName{2,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_2017\b01_2013-02-13_18-21-28';
% folderName{3,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\aveTreat_aveCatIV_empUpta_slowDST_2027\b01_2013-02-14_00-46-02';
% folderName{4,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\feb6_smokingUnstrat_vynnicky\base_60perc_private0p5\b2_2013-02-11_15-21-51';

% 
% folderName{1,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_1997\a03_2014-01-11_12-28-21';
% folderName{2,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2007\a03_2014-01-11_14-20-56';
% folderName{3,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2017\a03_2014-01-10_16-32-47';
% folderName{4,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2017\a03_2014-01-11_16-13-42';
% folderName{5,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2027\a03_2014-01-10_18-29-33';
% folderName{6,1} = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\Jan9\aveTreat_aveCatIV_empUpta_slowDST_2027\a03_2014-01-11_18-06-45';

%or you can do it automatically if it's easier
% masterFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\jan18\base60_0p55_rntcpMDR1mon\';
% foldersStruc = dir(masterFolder);
% is_dir = [foldersStruc(:).isdir]';
% foldnameArray_all = {foldersStruc(:).name}';
% foldName = foldnameArray_all(is_dir);
% foldName = folderName(3:end); % Exclude . and ..
% folderName = strcat(masterFolder,foldName{folds} );

for folds = 1:size(folderName,1)
    foldNameNow = folderName{folds};
    cd(foldNameNow);
    diagnosedPpl = dlmread('diagnosedPpl.csv', ',' ,1,0);
    diag2008 = sum(diagnosedPpl(1704:1715,:),1);
    incMat = dlmread('incidenceMatrix_forMakingWHOComparison.csv', ',' ,1,0);  
	
    popIndia2008 = 1190863679;
    alpha = 2.5;
%    inciRateInHunThou = 186;
    inciRateInHunThou = sum(incMat(1704:1715,1)./incMat(1704:1715,3))*100000;
    
    numMDR_CatI = diag2008(1,9);
    numPpl_CatI = diag2008(1,6);
    numRelapse = diag2008(1,11);
    numMDR_CatII = diag2008(1,10);
    num_CatII = diag2008(1,7);
    numMDR_relapsers = diag2008(1,12);
    
    m_n = numMDR_CatI/numPpl_CatI;
    inciIndianPop =inciRateInHunThou*(popIndia2008/100000);
    propMDR_CatII_excludingRelapsers =(numMDR_CatII-numMDR_relapsers)/(num_CatII-numRelapse);
    nonRelapseRetreat =num_CatII-numRelapse;
    newPlusRelapse =numRelapse+numPpl_CatI;
    MDRinAllTreatedCases =(numMDR_CatII+numMDR_CatI)/(numPpl_CatI+num_CatII);
    relapseInRetreat =numMDR_relapsers/num_CatII;
    relapseOverNewPlusRelapse =numRelapse/newPlusRelapse;
    
    I_m(folds) = ((m_n*(1-relapseOverNewPlusRelapse))+(alpha*m_n*relapseOverNewPlusRelapse))*inciIndianPop;
    I_o1(folds) = (propMDR_CatII_excludingRelapsers-m_n)*(nonRelapseRetreat/newPlusRelapse)*inciIndianPop;
    I_o2(folds) =(((MDRinAllTreatedCases-(alpha*relapseInRetreat*m_n))/(1-relapseInRetreat))-m_n)*(nonRelapseRetreat/newPlusRelapse)*inciIndianPop;
    sumI_o1(folds) = I_m(folds) + I_o1(folds);
    sumI_o2(folds) = I_m(folds) + I_o2(folds);
    percMDRNewCases(folds) = ((m_n*(1-relapseOverNewPlusRelapse))+(alpha*m_n*relapseOverNewPlusRelapse));
end


percMDRNewCases_means =  mean(percMDRNewCases,2);
[a,b,percMDRNewCases_CIs] = ttest(percMDRNewCases);
I_m_means =  mean(I_m,2);
[a,b,I_m_CIs] = ttest(I_m);
sumI_o1_means =  mean(sumI_o1,2);
[a,b,sumI_o1_CIs] = ttest(sumI_o1);
sumI_o2_means =  mean(sumI_o2,2);
[a,b,sumI_o2_CIs] = ttest(sumI_o2);
format longG
MDRoutVals = [percMDRNewCases_means,percMDRNewCases_CIs;...
I_m_means,I_m_CIs;...
sumI_o1_means,sumI_o1_CIs;...
sumI_o2_means,sumI_o2_CIs]

if makeGraphs == 1
    restricted = 0;
    noAndBigSeed = 0;
    
    subplot(1,4,1);
    percMDRinNewCases{1} = [...
        2008    2.3	1.8	2.8
        2010    2.1	1.5	2.7
        ];
    percMDRinNewCases{2}=[
        2008.2  1.4	1.0 1.8
        2010.2  1.6	1.2	2.0
        ];
    percMDRinNewCases{5} = [...
        2008.8 0.66620058   0.66620058  0.66620058
        2010.8    0.885889097 0.885889097   0.885889097
        ];
    
    percMDRinNewCases{6} = [...
        2008.8  7.237502657 7.237502657 7.237502657
        2010.8  7.427596136  7.427596136    7.427596136
        ];
    
    errorbar(percMDRinNewCases{1}(:,1),percMDRinNewCases{1}(:,2), percMDRinNewCases{1}(:,2)-percMDRinNewCases{1}(:,3), percMDRinNewCases{1}(:,4)-percMDRinNewCases{1}(:,2), 'bo');
    hold on;
    errorbar(percMDRinNewCases{2}(:,1),percMDRinNewCases{2}(:,2), percMDRinNewCases{2}(:,2)-percMDRinNewCases{2}(:,3), percMDRinNewCases{2}(:,4)-percMDRinNewCases{2}(:,2), 'go');
    if noAndBigSeed == 1
        hold on;
        errorbar(percMDRinNewCases{5}(:,1),percMDRinNewCases{5}(:,2), percMDRinNewCases{5}(:,2)-percMDRinNewCases{5}(:,3), percMDRinNewCases{5}(:,4)-percMDRinNewCases{5}(:,2), 'mo');
        hold on;
        errorbar(percMDRinNewCases{6}(:,1),percMDRinNewCases{6}(:,2), percMDRinNewCases{6}(:,2)-percMDRinNewCases{6}(:,3), percMDRinNewCases{6}(:,4)-percMDRinNewCases{6}(:,2), 'ko');
    end
    ylabel('Percent MDR among new TB cases');
    ylim([0 3.5]);
    xlim([2007 2009.9]);
    if noAndBigSeed == 1
        ylim([0 8]);
        xlim([2007 2011]);
    end
    set(gca,'XTick',[2008 2010]);
    
    %I_m
    subplot(1,4,2);
    percMDRinNewCases{1} = [...
        2008    55000	40000	74000
        2010    21000	15000	27000
        ];
    percMDRinNewCases{2}=[
        2008.2  24042	22872	26336
        2010.2  25206	20245	30166
        ];
    %equation 1
    percMDRinNewCases{3}=[
        2008.4  30796	26958	34633
        2010.4  35684	27000	44368
        ];
    %noSeed
    percMDRinNewCases{5} = [...
        2008.8 1475 1475    1475
        2010.8  19962   19962   19962
        ];
    %bigSeed
    percMDRinNewCases{6} = [...
        2008.8  160311  160311  160311
        2010.8  164522  164522  164522
        ];
    
    errorbar(percMDRinNewCases{1}(:,1),percMDRinNewCases{1}(:,2), percMDRinNewCases{1}(:,2)-percMDRinNewCases{1}(:,3), percMDRinNewCases{1}(:,4)-percMDRinNewCases{1}(:,2), 'bo');
    hold on;
    if restricted == 0
        errorbar(percMDRinNewCases{2}(:,1),percMDRinNewCases{2}(:,2), percMDRinNewCases{2}(:,2)-percMDRinNewCases{2}(:,3), percMDRinNewCases{2}(:,4)-percMDRinNewCases{2}(:,2), 'ro');
        hold on;
    end
    if noAndBigSeed == 1;
        errorbar(percMDRinNewCases{5}(:,1),percMDRinNewCases{5}(:,2), percMDRinNewCases{5}(:,2)-percMDRinNewCases{5}(:,3), percMDRinNewCases{5}(:,4)-percMDRinNewCases{5}(:,2), 'mo');
        hold on;
        errorbar(percMDRinNewCases{6}(:,1),percMDRinNewCases{6}(:,2), percMDRinNewCases{6}(:,2)-percMDRinNewCases{6}(:,3), percMDRinNewCases{6}(:,4)-percMDRinNewCases{6}(:,2), 'ko');
        hold on;
    end
    errorbar(percMDRinNewCases{3}(:,1),percMDRinNewCases{3}(:,2), percMDRinNewCases{3}(:,2)-percMDRinNewCases{3}(:,3), percMDRinNewCases{3}(:,4)-percMDRinNewCases{3}(:,2), 'go');
    ylabel({'Number of MDR-TB among incident new and relapse TB cases'});
    xlim([2007 2009.9]);
    ylim([0 80000]);
    if restricted == 0 || noAndBigSeed == 1
        ylim([0 170000]);
        xlim([2007 2011]);
    end
    set(gca,'XTick',[2008 2010]);
    
    subplot(1,4,3);
    percMDRinNewCases{1} = [...
        2008    43000	33000	56000
        2010    43000	39000	48000
        ];
    percMDRinNewCases{2}=[
        2008.2  55496	50431	65422
        2010.2  54827	42951	66704
        ];
    %equation 1
    percMDRinNewCases{3}=[
        2008.4  263363	196502  330224
        2010.4  266703	233188  300219
        ];
    %equation 2
    percMDRinNewCases{4}=[
        2008.6  66796	41325   92267
        2010.6  65444	51390   79497
        ];
    percMDRinNewCases{5} = [...
        2008.8 244787   244787  244787
        2010.8 253786   244787  244787
        ];
    percMDRinNewCases{6} = [...
        2008.8  412788  412788  412788
        2010.8  387341  387341  387341
        ];
    percMDRinNewCases{7} = [...
        2008.8 58262    58262   58262
        2010.8 	60721   60721   60721
        ];
    
    percMDRinNewCases{8} = [...
        2008.8  140957  140957  140957
        2010.8  123374  140957 140957
        ];
    
    
    errorbar(percMDRinNewCases{1}(:,1),percMDRinNewCases{1}(:,2), percMDRinNewCases{1}(:,2)-percMDRinNewCases{1}(:,3), percMDRinNewCases{1}(:,4)-percMDRinNewCases{1}(:,2), 'bo');
    hold on;
    if restricted == 0
        errorbar(percMDRinNewCases{2}(:,1),percMDRinNewCases{2}(:,2), percMDRinNewCases{2}(:,2)-percMDRinNewCases{2}(:,3), percMDRinNewCases{2}(:,4)-percMDRinNewCases{2}(:,2), 'ro');
        hold on;
        errorbar(percMDRinNewCases{3}(:,1),percMDRinNewCases{3}(:,2), percMDRinNewCases{3}(:,2)-percMDRinNewCases{3}(:,3), percMDRinNewCases{3}(:,4)-percMDRinNewCases{3}(:,2), 'gx');
        hold on;
    end
    if noAndBigSeed == 1;
        errorbar(percMDRinNewCases{5}(:,1),percMDRinNewCases{5}(:,2), percMDRinNewCases{5}(:,2)-percMDRinNewCases{5}(:,3), percMDRinNewCases{5}(:,4)-percMDRinNewCases{5}(:,2), 'mx');
        hold on;
        errorbar(percMDRinNewCases{6}(:,1),percMDRinNewCases{6}(:,2), percMDRinNewCases{6}(:,2)-percMDRinNewCases{6}(:,3), percMDRinNewCases{6}(:,4)-percMDRinNewCases{6}(:,2), 'kx');
        hold on;
        errorbar(percMDRinNewCases{7}(:,1),percMDRinNewCases{7}(:,2), percMDRinNewCases{7}(:,2)-percMDRinNewCases{7}(:,3), percMDRinNewCases{7}(:,4)-percMDRinNewCases{7}(:,2), 'mo');
        hold on;
        errorbar(percMDRinNewCases{8}(:,1),percMDRinNewCases{8}(:,2), percMDRinNewCases{8}(:,2)-percMDRinNewCases{8}(:,3), percMDRinNewCases{8}(:,4)-percMDRinNewCases{8}(:,2), 'ko');
        hold on;
    end
    errorbar(percMDRinNewCases{4}(:,1),percMDRinNewCases{4}(:,2), percMDRinNewCases{4}(:,2)-percMDRinNewCases{4}(:,3), percMDRinNewCases{4}(:,4)-percMDRinNewCases{4}(:,2), 'go');
    ylabel({'Number of incident acquired MDR-TB cases'});
    xlim([2007 2009.9]);
    ylim([0 100000]);
    if restricted == 0 || noAndBigSeed == 1
        xlim([2007 2011]);
        ylim([0 500000]);
    end
    set(gca,'XTick',[2008 2010]);
    
    
    subplot(1,4,4);
    percMDRinNewCases{1} = [...
        2008    99000	79000	120000
        2010    64000   54000   75000
        ];
    percMDRinNewCases{2}=[
        2008.2  79538	74265	89873
        2010.2  80033	66446	93620
        ];
    %equation 1
    percMDRinNewCases{3}=[
        2008.4  294159	227241  361077
        2010.4  302387	274214  330561
        ];
    %equation 2
    percMDRinNewCases{4}=[
        2008.6  97591	71969   123213
        2010.6  101128	92242   110014
        ];
    percMDRinNewCases{5} = [...
        2008.8 259544   259544  259544
        2010.8 273748   273748  273748
        ];
    
    percMDRinNewCases{6} = [...
        2008.8  573099   573099   573099
        2010.8  551863  551863    551863
        ];
    %total_2
    percMDRinNewCases{7} = [...
        2008.8 73018    73018    73018
        2010.8 80682    80682    80682
        ];
    
    percMDRinNewCases{8} = [...
        2008.8 301268    301268   301268
        2010.8  287896   287896   287896
        ];
    
    errorbar(percMDRinNewCases{1}(:,1),percMDRinNewCases{1}(:,2), percMDRinNewCases{1}(:,2)-percMDRinNewCases{1}(:,3), percMDRinNewCases{1}(:,4)-percMDRinNewCases{1}(:,2), 'bo');
    hold on;
    if restricted == 0
        errorbar(percMDRinNewCases{2}(:,1),percMDRinNewCases{2}(:,2), percMDRinNewCases{2}(:,2)-percMDRinNewCases{2}(:,3), percMDRinNewCases{2}(:,4)-percMDRinNewCases{2}(:,2), 'ro');
        hold on;
        errorbar(percMDRinNewCases{3}(:,1),percMDRinNewCases{3}(:,2), percMDRinNewCases{3}(:,2)-percMDRinNewCases{3}(:,3), percMDRinNewCases{3}(:,4)-percMDRinNewCases{3}(:,2), 'gx');
        hold on;
    end
    if noAndBigSeed == 1;
        errorbar(percMDRinNewCases{5}(:,1),percMDRinNewCases{5}(:,2), percMDRinNewCases{5}(:,2)-percMDRinNewCases{5}(:,3), percMDRinNewCases{5}(:,4)-percMDRinNewCases{5}(:,2), 'mx');
        hold on;
        errorbar(percMDRinNewCases{6}(:,1),percMDRinNewCases{6}(:,2), percMDRinNewCases{6}(:,2)-percMDRinNewCases{6}(:,3), percMDRinNewCases{6}(:,4)-percMDRinNewCases{6}(:,2), 'kx');
        hold on;
        errorbar(percMDRinNewCases{7}(:,1),percMDRinNewCases{7}(:,2), percMDRinNewCases{7}(:,2)-percMDRinNewCases{7}(:,3), percMDRinNewCases{7}(:,4)-percMDRinNewCases{7}(:,2), 'mo');
        hold on;
        errorbar(percMDRinNewCases{8}(:,1),percMDRinNewCases{8}(:,2), percMDRinNewCases{8}(:,2)-percMDRinNewCases{8}(:,3), percMDRinNewCases{8}(:,4)-percMDRinNewCases{8}(:,2), 'ko');
        hold on;
    end
    rows = [1 2];
    errorbar(percMDRinNewCases{4}(rows,1),percMDRinNewCases{4}(rows,2), percMDRinNewCases{4}(rows,2)-percMDRinNewCases{4}(rows,3), percMDRinNewCases{4}(rows,4)-percMDRinNewCases{4}(rows,2), 'go');
    ylabel({'Number of MDR-TB among incident total TB cases' });
    xlim([2007 2009.9]);
    ylim([0 150000]);
    if restricted == 0 || noAndBigSeed == 1
        xlim([2007 2011]);
        ylim([0 500000]);
    end
    set(gca,'XTick',[2008 2010]);
    if restricted == 0
        legend('WHO','Simulation: directly observed','Simulation: Calculated with WHO eqnA','Simulation: Calculated with WHO eqnB' ,'Location', 'NorthEast');
    end
    if noAndBigSeed ==1
        legend('WHO', 'Simulation, No MDR Seed: Calculated','Simulation, Max MDR Seed: Calculated' ,'Simulation: Calculated with WHO eqn1','Location', 'NorthEast');
    end
    plotResolution = '-r100';
    outFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\jan14\paperFigures\mdr_whoComp';
    print('-dpng',plotResolution,[outFold '/'  'who_MDRcomparison']);
    
    
end