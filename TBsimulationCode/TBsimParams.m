%TBparams
%TBparams.percAgeBrac
%TBparams.meanAge
%TBparams.ageBracPercMale
%TBparams.ruralTBprev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Year to month num Translation%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TBparams.yearToMon = [[1866:1:2045]', [1:12:2149]'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BASELINE AGE AND GENDER%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TBparams.cumpercAgeBrac = [0;...
0.073165976 ;...
0.141385932 ;...
0.207752959 ;...
0.273420886 ;...
0.338619453 ;...
0.403128065 ;...
0.466729419 ;...
0.52930462  ;...
0.590794769 ;...
0.650987686 ;...
0.709607168 ;...
0.7660485   ;...
0.819340767 ;...
0.867916756 ;...
0.910525328 ;...
0.944774687 ;...
0.970026826 ;...
0.985756414 ;...
0.994406827 ;...
0.99827621  ;...
0.999644048 ;...
1   ;...
];
TBparams.cumpercAgeBracLimits = [...
0   ;... 
1   ;... 
5   ;... 
10  ;... 
15  ;... 
20  ;... 
25  ;... 
30  ;... 
35  ;... 
40  ;... 
45  ;... 
50  ;... 
55  ;... 
60  ;... 
65  ;... 
70  ;... 
75  ;... 
80  ;... 
85  ;... 
90  ;... 
95  ;... 
99  ;... 
100 ;...
];
TBparams.percMaleAgeBrac = [...
0   ,   1   ;... 
1   ,   4   ;... 
5   ,   9   ;... 
10  ,   14  ;... 
15  ,   19  ;... 
20  ,   24  ;... 
25  ,   29  ;... 
30  ,   34  ;... 
35  ,   39  ;... 
40  ,   44  ;... 
45  ,   49  ;... 
50  ,   54  ;... 
55  ,   59  ;... 
60  ,   64  ;... 
65  ,   69  ;... 
70  ,   74  ;... 
75  ,   79  ;... 
80  ,   84  ;... 
85  ,   89  ;... 
90  ,   94  ;... 
95  ,   99  ;... 
100 ,   Inf ;... 
];
TBparams.percMale = [0.5    ;...
0.500380738 ;...
0.503282529 ;...
0.503916347 ;...
0.503981012 ;...
0.50468143  ;...
0.505375169 ;...
0.505401929 ;...
0.504432308 ;...
0.502248721 ;...
0.499335359 ;...
0.494973523 ;...
0.489308392 ;...
0.480972715 ;...
0.470228641 ;...
0.45599812  ;...
0.442609414 ;...
0.421099146 ;...
0.400490569 ;...
0.379880874 ;...
0.362396363 ;...
0.348406989 ;...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%BASELINE LATENT TB PREVALENCE%%%%%%%%%%%%%
TBparams.latentTBageBrac = [...
0, 4;...
5, 9;...
10, 14;...
15, 19;...
20, 24;...
25, 29;...
30, 34;...
35, 39;...
40, 44;...
45, 49;...
50, 54;...
55, 59;...
60, 64;...
65, 69;...
70, 74;...
75, 79;...
80, 84;...
85, 89;...
90, 94;...
95, Inf;...
];
%this is the original latent TB distribution, but makes latent too high in
%1996 so need to scale down.
% TBparams.latentTBprev=[...
% .11;...
% .22;...
% .33;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% .4399494469;...
% ];
TBparams.latentTBprev=[...
.04;...
.08;...
.08;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
.0899494469;...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%BASELINE TB PREVALENCE%%%%%%%%%%%%%
%columns: male, female
%rows:
TBparams.ruralTBprevAgeBrac=[15,24;...
25,34;...
35,44;...
45,54;...
55,64;...
65,Inf...
];
TBparams.ruralTBprev = [...
0.00249,    0.00047;...
0.00693,    0.00223;...
0.01527,    0.00303;...
0.02681,    0.00279;...
0.04313,    0.00681;...
0.04534,    0.00605;...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prevalence of TB higher for smokers: this is only used in initial cohort
TBparams.TBprevSmokingMultiplier = 3.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%SMOKING PREVALENCE%%%%%%%%%%%%%
TBparams.smokingPrevAgeBrac = [...  
0   ,   17  ;...
18  ,   24  ;...
25  ,   29  ;...
30  ,   34  ;...
35  ,   39  ;...
40  ,   44  ;...
45  ,   49  ;...
50  ,   54  ;...
55  ,   59  ;...
60  ,   64  ;...
65  ,   Inf ;...
];
TBparams.smokingPrev =[...
    0   ,   0    ;...
0.144027,0.020664;...
0.302404,0.027755;...
0.381743,0.038614;...
0.415032,0.049292;...
0.430601,0.036372;...
0.424791,0.05139;...
0.413023,0.070911;...
0.406085,0.054887;...
0.381571,0.059954;...
0.300339,0.047119;...
];
%%%%%%%%%%%%%%%SMOKING CHURN%%%%%%%%%%%%%%%%%
%old numbers not using since smoking does not look right
% TBparams.smokingChurnAgeBrac = [5:1:100;6:1:101]';
% TBparams.smoker_quit_ratio = 0.464159;
% TBparams.nonsmoker_start_ratio = 0.736806;
%jeremy numbers after 1pr24
TBparams.smokingChurnAgeBrac = [5:1:100;6:1:101]';
TBparams.smoker_quit_ratio=0.577981;
TBparams.nonsmoker_start_ratio=0.468645;
[TBparams.oneMonQuittingUrbanProb, TBparams.oneMonQuittingRuralProb, TBparams.oneMonStartingUrbanProb, TBparams.oneMonStartingRuralProb] = smokingChurnMaker;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Population growth over 10 years
TBparams.popGrowthPerc = 0.1165;
TBparams.hisPopGrowthYrs = [...
1995
2000
2005
2010
2015
2020
2025
2030
2035
2040
2045
2050    
2055
2060
2065
2070
2075
2080
2085
2090
2095
2100
];
%gives historic population growth trend (rather than constant growth rate)
TBparams.hisPopGrowthPerc = [...
0.020760485
0.01854086
0.016347882
0.014836458
0.013654425
0.012029772
0.010389867
0.008845217
0.007393589
0.005978851
0.0046084
0.003302936
0.002095144
0.000963425
-8.97571E-05
-0.001047986
-0.00187238
-0.002570133
-0.003132307
-0.003566989
-0.003897648
-0.004117711
];
%gives a growth rate such that will always give larger number if compound
%at historic growth rates over duration of simulation (using annual compounding
%for both historical and constant rates).  This is used to make the total
%number of possible people (length of stateMat).  Determined through guess
%and check using the historic data.
TBparams.constantHisPopGrowth = 0.018183;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%historic population in India
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Base mortality
TBparams.fixedDeathYr = 0;
TBparams.baseMortAgeBrac = [...
0   ,   0   ;...
1   ,   4   ;...
5   ,   9   ;...
10  ,   14  ;...
15  ,   19  ;...
20  ,   24  ;...
25  ,   29  ;...
30  ,   34  ;...
35  ,   39  ;...
40  ,   44  ;...
45  ,   49  ;...
50  ,   54  ;...
55  ,   59  ;...
60  ,   64  ;...
65  ,   69  ;...
70  ,   74  ;...
75  ,   79  ;...
80  ,   84  ;...
85  ,   89  ;...
90  ,   94  ;...
95  ,   98  ;...
99  ,   Inf ];

datayears = [1990; 2000; 2009];


mortMaleOverTime = [...
0.007313128	0.005831265	0.004287449
0.000645625	0.000454064	0.000259966
0.000233306	0.000155821	0.00011416
0.000137491	0.000117493	8.08301E-05
0.000166653	0.000154155	0.000131658
0.000230807	0.000213311	0.00018415
0.000255801	0.000269964	0.00021831
0.00030412	0.000323281	0.000278295
0.000385759	0.000427409	0.000379095
0.000531525	0.000538188	0.000499042
0.000865459	0.000777198	0.000690595
0.001269194	0.001147674	0.00100283
0.001985526	0.001827495	0.001413999
0.003024583	0.00255257	0.002412086
0.004503993	0.004122313	0.003634214
0.006916801	0.005508937	0.005701186
0.0096333	0.008489591	0.008081336
0.013190565	0.010384868	0.011574165
0.018051746	0.013383796	0.01643921
0.024686837	0.018170391	0.023149508
0.033728086	0.025971778	0.032309847
1 1 1
];

mortFemaleOverTime=[...
0.007477735	0.005960499	0.004382037
0.000965367	0.000697257	0.000425743
0.000290791	0.00019748	0.00011916
0.000148322	0.000121659	8.49964E-05
0.000242471	0.000200813	0.000131658
0.000309952	0.000259133	0.000176651
0.000290791	0.00027163	0.000164986
0.000274962	0.000259133	0.000186649
0.000329112	0.00028246	0.000219976
0.000401586	0.000344107	0.000280794
0.000572336	0.000486548	0.000382427
0.000870454	0.000770536	0.000561509
0.001352418	0.001273355	0.00096953
0.002248302	0.00183914	0.001750965
0.003635874	0.003183256	0.002691372
0.005884286	0.004633399	0.004661601
0.007918482	0.007120363	0.006536874
0.012395864	0.009100834	0.010031842
0.018586762	0.012170335	0.015037458
0.026689859	0.017015244	0.022012444
0.036704442	0.024864003	0.031457098
1 1 1
];

mortMale = interp1(datayears, mortMaleOverTime', [1990:1:2009]);
mortFemale = interp1(datayears, mortFemaleOverTime', [1990:1:2009]);
for i = 1 : ( max(datayears) - min(datayears) + 1)
    TBparams.mortNonsmoking{i} = [mortMale(i,:)', mortFemale(i,:)'];
end


%took out the smoking mix, now it's just the who life tables
% 
% mortMaleSmokingOverTime = [...
% 0.007313128 0.005831265 0.004287449
% 0.000645625 0.000454064 0.000259966
% 0.000233306 0.000155821 0.00011416
% 0.000137491 0.000117493 8.08301E-05
% 0.000166653 0.000154155 0.000131658
% 0.000269213 0.000248806 0.000214792
% 0.000289454 0.000305481 0.000247031
% 0.000355906 0.00038562  0.000333578
% 0.000508908 0.00057271  0.000509989
% 0.000695301 0.000714684 0.000665238
% 0.001132751 0.001032905 0.000921383
% 0.00167174  0.001535609 0.001347174
% 0.002532309 0.002369802 0.001841369
% 0.003901856 0.003351088 0.0031808
% 0.00583019  0.005449971 0.004830619
% 0.00945267  0.007528648 0.007791381
% 0.0096333   0.008489591 0.008081336
% 0.013190565 0.010384868 0.011574165
% 0.018051746 0.013383796 0.01643921
% 0.024686837 0.018170391 0.023149508
% 0.033728086 0.025971778 0.032309847
% 1   1   1
% ];
% mortFemaleSmokingOverTime=[...
% 0.007477735 0.005960499 0.004382037
% 0.000965367 0.000697257 0.000425743
% 0.000290791 0.00019748  0.00011916
% 0.000148322 0.000121659 8.49964E-05
% 0.000242471 0.000200813 0.000131658
% 0.000460173 0.000384725 0.000262267
% 0.000430216 0.000401868 0.000244092
% 0.000370422 0.000358287 0.000259085
% 0.000564643 0.000496911 0.000388458
% 0.00069575  0.000611674 0.000501073
% 0.000975154 0.00084997  0.000670606
% 0.001461918 0.001325702 0.000969601
% 0.002144436 0.002070207 0.001582224
% 0.00355382  0.002980056 0.002847839
% 0.005486522 0.004927354 0.004182038
% 0.007555191 0.005949102 0.005985312
% 0.007918482 0.007120363 0.006536874
% 0.012395864 0.009100834 0.010031842
% 0.018586762 0.012170335 0.015037458
% 0.026689859 0.017015244 0.022012444
% 0.036704442 0.024864003 0.031457098
% 1   1   1
% ];
% mortMaleNonsmokingOverTime=[...
% 0.007313128 0.005831265 0.004287449
% 0.000645625 0.000454064 0.000259966
% 0.000233306 0.000155821 0.00011416
% 0.000137491 0.000117493 8.08301E-05
% 0.000166653 0.000154155 0.000131658
% 0.000224344 0.000207338 0.000178994
% 0.000241212 0.000254567 0.000205859
% 0.000272146 0.00028479  0.00024416
% 0.000298385 0.000324318 0.000286226
% 0.000407672 0.000404715 0.000373358
% 0.000668064 0.000588358 0.000520158
% 0.000985945 0.000874706 0.000760534
% 0.001611668 0.001456696 0.001121789
% 0.002483306 0.002059884 0.00193779
% 0.003934705 0.003552398 0.00312064
% 0.005907919 0.004705405 0.004869613
% 0.0096333   0.008489591 0.008081336
% 0.013190565 0.010384868 0.011574165
% 0.018051746 0.013383796 0.01643921
% 0.024686837 0.018170391 0.023149508
% 0.033728086 0.025971778 0.032309847
% 1   1   1
% ];
% mortFemaleNonsmokingOverTime=[...
% 0.007477735 0.005960499 0.004382037
% 0.000965367 0.000697257 0.000425743
% 0.000290791 0.00019748  0.00011916
% 0.000148322 0.000121659 8.49964E-05
% 0.000242471 0.000200813 0.000131658
% 0.000306782 0.000256483 0.000174845
% 0.000286811 0.000267912 0.000162728
% 0.000271128 0.000255151 0.00018374
% 0.000316901 0.000271341 0.00021124
% 0.000390483 0.000334008 0.00027248
% 0.000550514 0.00046686  0.000366815
% 0.000825312 0.000728164 0.000530362
% 0.001306422 0.001227078 0.000933948
% 0.002165039 0.001766375 0.001681009
% 0.003544361 0.003097012 0.00261766
% 0.005811686 0.004576232 0.004604086
% 0.007918482 0.007120363 0.006536874
% 0.012395864 0.009100834 0.010031842
% 0.018586762 0.012170335 0.015037458
% 0.026689859 0.017015244 0.022012444
% 0.036704442 0.024864003 0.031457098
% 1   1   1
% % ];
% mortMaleSmoking = interp1(datayears, mortMaleSmokingOverTime', [1990:1:2009]);
% mortFemaleSmoking = interp1(datayears, mortFemaleSmokingOverTime', [1990:1:2009]);
% mortMaleNonsmoking = interp1(datayears, mortMaleNonsmokingOverTime', [1990:1:2009]);
% mortFemaleNonsmoking = interp1(datayears, mortFemaleNonsmokingOverTime', [1990:1:2009]);
% for i = 1 : ( max(datayears) - min(datayears) + 1)
%     TBparams.mortSmoking{i} = [mortMaleSmoking(i,:)', mortFemaleSmoking(i,:)'];
%     TBparams.mortNonsmoking{i} = [mortMaleNonsmoking(i,:)', mortFemaleNonsmoking(i,:)'];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fraction  of population in urban in 2001
%male   %female
TBparams.urbanFrac = [...
0.282594338,    0.27276306...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mortality prob from non-treated TB;
untreatedTBmortProbNonSmoker = 0.02469;  %this is actually a rate
%smokingTBmortRiskRatio = 1.475;  %not doing smoking anymore
%untreatedTBmortProbSmoker = untreatedTBmortProbNonSmoker * smokingTBmortRiskRatio;  %not doing smoking anymore
warning off
untreatedTBmortRateNonSmok = -log(1-untreatedTBmortProbNonSmoker);
%untreatedTBmortRateSmok = -log(1-untreatedTBmortProbSmoker);
warning on
for i = 1: ( max(datayears) - min(datayears) + 1)
    warning off
    %baseMortSmokingRate{i} = -log(1-TBparams.mortSmoking{i});
    baseMortNonsmokingRate{i} = -log(1-TBparams.mortNonsmoking{i});
    warning on
    %TBparams.nonDOTsTBmortSmoking{i} = 1-exp(-(baseMortSmokingRate{i} + untreatedTBmortRateSmok));
    TBparams.nonDOTsTBmortNonSmok{i} = 1-exp(-(baseMortNonsmokingRate{i} + untreatedTBmortRateNonSmok));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probability of self cure should be 0 in base case
TBparams.probSelfCure = 0.0064;  %for 20 percent over 3 years.  Base case used to be value = 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%want to seed population with treated TB ppl
%Pr(retreatment | TB ) = Pr(retreatment,TB)/Pr(TB) 
%                   = Pr(retreat,TB | tested)*Pr(tested) / Pr(TB) 
%                   = Pr(retreatment | tested)Pr(tested)/Pr(TB)
%                   = 0.011985761 for 1 quarter (2009Q1).  Ie, must be
%                   higher than this. (this is lower bound)
TBparams.alreadyTreatedB4 = 0.03;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%case detection probability for sputum smear process;
%from mortality page;
% WHOcaseDetectionProb =  0.0279357;  %this was validated by RNTCP 2009Q1 numbers
% SSsensit = 0.6; %SS test has 50% to 70% sensitivity; 
% TBparams.SScaseDetectionProb = SSsensit*WHOcaseDetectionRate;
% TBparams.ssforEveryCorrectTBcaseNumIncorrect=((1-SSsensit)/SSsensit);
% TBparams.rapidDiagCaseDetectionProb =  rapidDiagSens*WHOcaseDetectionRate;
% TBparams.rpdDiagforEveryCorrectTBcaseNumIncorrect=((1-rapidDiagSens)/rapidDiagSens);
%testing and diagnosis parameters (from treatment and test page)
%pr(tested | TB ) = Pr(TB|tested)*Pr(tested) / Pr(TB) 
%                = 1/current sensitivty * (totPopProbTested / current prev)
%                = 1/current sensitivity * testingCoefficient
%
%  Used in simulation: TestSensivity * pr(tested | TB)
% 
% testingCoefVec = [...
% 0   ;...
% 0.001128263 ;...
% 0.002256526 ;...
% 0.003384789 ;...
% 0.004513053 ;...
% 0.005641316 ;...
% 0.006769579 ;...
% 0.007897842 ;...
% 0.009026105 ;...
% 0.010154368 ;...
% 0.011282632 ;...
% 0.012410895 ;...
% 0.013539158 ;...
% 0.014667421 ;...
% 0.015795684 ;...
% 0.016923947 ;...
% 0.018052211 ;...
% 0.019180474 ;...
% 0.020308737 ;...
% 0.021437    ;...
% 0.021437    ;...
% ];
% 
% TBparams.prTestedGivenTBVec = (1/TBparams.SSsensit)*testingCoefVec;  
TBparams.SSsensit = 0.6; %SS test has 50% to 70% sensitivity; base is 60% 
TBparams.SSspec = 1;% 0.9840;  %changed this to see if it makes the difference
% TBparams.rapidDiagSens = 0.863;  %this is not used
   testingCoefficient = 1;  %this is the prob of being successfully entered into treatment given have actTB in 2010Q3
    TBparams.averageKnowledgeOfTB = aveUptake_cal; %this is the proportion of sick people who know about treatment, and it was calibrated
    TBparams.prTestedGivenTB = testingCoefficient*(1/TBparams.SSsensit)*(1/TBparams.averageKnowledgeOfTB); %this is the average "willingness to be tested"
    TBparams.priorTreatBoostFactor = cat2uptake_cal;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %selection: some people never do RNCTP
    TBparams.probSeekRNTCPtreatment = 1; %changed from 0.60;
    TBparams.RNTCPqueueTime = 1;  %not used, is hard coded as 1
    
%knowledge of tb for treatment uptake
    %knowledge of tb for treatment uptake
    TBparams.uptakeAgeBracs = [...
        0	4
        5	9
        10	14
        15	19
        20	24
        25	29
        30	34
        35	39
        40	44
        45	49
        50	54
        55	59
        60	64
        65	69
        70	74
        75	79
        80	84
        85	Inf
        ];
    %scaling factors to calibrate treatment uptake
    %men
    
    %%% old setting
    %TBparams.uptakeUrbanKnowledge = [...
    %    0.0282572	0.0265703
    %    0.03421145	0.034611175
    %    0.0391262	0.0411383
    %    0.04314695	0.046267925
    %    0.0464192	0.0501163
    %    0.04908845	0.052799675
    %    0.0513002	0.0544343
    %    0.05319995	0.055136425
    %    0.0549332	0.0550223
    %    0.05664545	0.054208175
    %    0.0584822	0.0528103
    %    0.06058895	0.050944925
    %    0.0631112	0.0487283
    %    0.06619445	0.046276675
    %    0.0699842	0.0437063
    %    0.07462595	0.041133425
    %    0.0802652	0.0386743
    %    0.08704745	0.036445175
    %    ];
    %%%
    
    %%% L setting (original setting)
    TBparams.uptakeUrbanKnowledge = [...
        0.0282572	0.0265703
        0.0282572	0.0265703
        0.03421145	0.034611175
        0.0391262	0.0411383
        0.04314695	0.046267925
        0.0464192	0.0501163
        0.04908845	0.052799675
        0.13294343	0.12393383 % 30 - 34
        0.1225361	0.12221708 % 35 - 39
        0.08431062	0.08530207 % 40 - 44
        0.08464646	0.08522047 % 45 - 49
        0.07465641	0.08945654 % 50 - 54
        0.09882135	0.11919841 % 55 - 59
        0.0631112	0.0487283
        0.06619445	0.046276675
        0.0699842	0.0437063
        0.07462595	0.041133425
        0.0802652	0.0386743
        0.08704745	0.036445175
        ];
    
    %%% testing setting sfw (test the effectiveness of this parameter)
    TBparams.testUptakeUrbanKnowledge = [...
        0.0282572	0.0265703
        0.03421145	0.034611175
        0.0391262	0.0411383
        0.04314695	0.046267925
        0.0464192	0.0501163
        0.04908845	0.052799675
        0.26	0.24 % 30 - 34
        0.24	0.244 % 35 - 39
        0.17	0.17 % 40 - 44
        0.17	0.17 % 45 - 49
        0.15	0.179 % 50 - 54
        0.196	0.24 % 55 - 59
        0.0631112	0.0487283
        0.06619445	0.046276675
        0.0699842	0.0437063
        0.07462595	0.041133425
        0.0802652	0.0386743
        0.08704745	0.036445175
        ];
    
    %%% new one for Bryan's code
    %%% new setting sfw (optimal for 25 years, K = 0.3222222)
    TBparams.afterUptakeUrbanKnowledge = [...
        0.0282572	0.0265703
        0.03421145	0.034611175
        0.0391262	0.0411383
        0.04314695	0.046267925
        0.0464192	0.0501163
        0.04908845	0.052799675
        0.13294343	0.12393383 % 30 - 34
        0.1225361	0.12221708 % 35 - 39
        0.08431062	0.08530207 % 40 - 44
        0.08464646	0.08522047 % 45 - 49
        0.07465641	0.13945654 % 50 - 54
        0.14882135	0.11919841 % 55 - 59
        0.0631112	0.0487283
        0.06619445	0.046276675
        0.0699842	0.0437063
        0.07462595	0.041133425
        0.0802652	0.0386743
        0.08704745	0.036445175
        ];
    
    
    TBparams.uptakeRuralKnowledge = TBparams.uptakeUrbanKnowledge;
    TBparams.catIIuptakeKnowledge = TBparams.priorTreatBoostFactor*TBparams.uptakeUrbanKnowledge;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%private clinic parameters
%TBparams.private_SSsensit = 0.5*TBparams.SSsensit is specified above
TBparams.SSsensitPriv = TBparams.SSsensit;  %0.6;  %same as base case public SSsensit
TBparams.SSspecPrivate = TBparams.SSspec;  %1;  %same as the base case SSspec
%TBparams.seekPrivate = 0.57;   %this is not one since we know 50%-70% of treatment seekers seek in private, so if only 60% of act seek rntcp, less than all of eligible private seekers enter private
TBparams.seekPrivate = 0;   %this is not one since we know 50%-70% of treatment seekers seek in private, so if only 60% of act seek rntcp, less than all of eligible private seekers enter private
%TBparams.maxPrivateTrtCount = 7;  %the max number of private clinics you seek
TBparams.maxPrivateTrtCount = 0;  %the max number of private clinics you seek
TBparams.privateClinicInsteadOfRNTCP = 0;  %this is for the scenario where RNTCP is replaced by private clinics.  1 to turn on scenario. 
TBparams.privateReferToRNTCP = 0; %public-private mix.  Prob of private referring patient to RNTCP.
TBparams.PPMeffectiveness = 1;  % = 1 means that PPM refers everyone it touches

TBparams.privateClinicCure = 0.0211;  %prob of cure in private clinic, for scenario for when private clinic can cure.  default 0.  use 0.0211?
TBparams.privateMDRprob = 0.02277884 * 0.242; %default is 0.02277884 * 0.242 which is the max default prob of all age/sex grps * TBparams.prbMDRafterDefault

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DOTS Plus Scale up!
%scales up dots plus coverage (yearly)
%I don't seem to use this anywhere
% TBparams.dotsPlusScaleUp = [...
% 0   ;...
% 0.180207447 ;...
% 0.360414893 ;...
% 0.54062234  ;...
% 0.720829786 ;...
% 0.901037233 ;...
% 1   ;...
% ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cat I probs;
TBparams.TBCatIageBrac = [
 0  19
20  29
30  39
40  49
50  59
60  69
70  Inf
];
TBparams.prbCatIrelapse = 0.197;
%NOT DOING THIS ANYMORE
%%%%%%%The "time zero" treatment quality, for the treatment quality ramp up%%%%%%%%%%%%%%%%%
%these are the "worst" case scenarios, since I don't have good data on
%historic treatment.  These are linearly ramped up over the RNTCP ramp up
%period to achieve the "average" (or "best" or "worst") treatment quality
%by 2007
% % %male   %female
% TBparams.timeZeroCatIdeath = [...
%     0.04645828  ,   0.04645828  ;...
%     0.04645828  ,   0.04645828  ;...
%     ];
% %given not dead
% TBparams.timeZeroCatIdefault = [...
%     0.020613709   ,   0.008808119 ;...
%     0.08584862    ,   0.036682621 ;...
%     ];
% %given have an outcome (ie, time to have outcome, not dead or defaulted)
% TBparams.timeZeroCatIsuccess = [...
%     0.90  ,   0.90    ;...
%     0.90  ,   0.90    ;...
%     ];
%%%%%%End the "time zero" treatment quality
%@@@@@@@@@@@using the AVERAGE parameters
%male   %female
     TBparams.CatIdeath = repmat([ 0.010101701 ,   0.010101701 ],size(TBparams.TBCatIageBrac,1),1);
     %given not dead
     TBparams.CatIdefault = [...
         0.0229393	0.0228821
         0.0222473	0.0191721
         0.0219013	0.0173171
         0.0215553	0.0154621
         0.0212093	0.0136071
         0.0208633	0.0117521
         0.0205173	0.0098971
         ];
     %given have an outcome (ie, time to have outcome, not dead or defaulted)
     TBparams.CatIsuccess = repmat([ 0.98   ,   0.98 ],size(TBparams.TBCatIageBrac,1),1);
   %WHO thinks these are closer to 88%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cat II death prob;
TBparams.TBCatIIageBrac = [
 0  19
20  29
30  39
40  49
50  59
60  69
70  Inf
];

TBparams.prbCatIIrelapse = 0.197;
%NOT USED ANYMORE
%%%%%%%The "time zero" treatment quality, for the treatment quality ramp up%%%%%%%%%%%%%%%%%
%these are the "worst" case scenarios, since I don't have good data on
%historic treatment.  These are linearly ramped up over the RNTCP ramp up
%period to achieve the "average" (or "best" or "worst") treatment quality
% %by 2007
% %male   %female
% TBparams.timeZeroCatIIdeath = [...
%     0.04135318  ,   0.04135318  ;...
%     0.04135318  ,   0.04135318  ;...
%     ];
% %given not deads
% TBparams.timeZeroCatIIdefault = [...
%     0.102510253   ,   0.048234838 ;...
%     0.1473321 ,   0.069325162 ;...
%     ];
% %given have an outcome (ie, time to have outcome, not dead or defaulted)
% TBparams.timeZeroCatIIsuccess = [...
%     0.83  ,   0.83    ;...
%     0.83  ,   0.83    ;...
%     ];
%%%%%%End the "time zero" treatment quality
%@@@@@@@@@@@using the AVERAGE parameters
     %CAT II
     TBparams.CatIIdeath = repmat([ 0.026003663 ,   0.026003663],size(TBparams.TBCatIageBrac,1),1);
     TBparams.CatIIdefault = [...
         0.0558639	0.0557247
         0.0541799	0.0466887
         0.0533379	0.0421707
         0.0524959	0.0376527
         0.0516539	0.0331347
         0.0508119	0.0286167
         0.0499699	0.0240987
         ];
     TBparams.CatIIsuccess = repmat([0.94   ,   0.94],size(TBparams.TBCatIageBrac,1),1);
TBparams.testPos4MonCatII = 0.57;  %from article 6 (see treatment tab, search CATIV TESTING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cat IV death prob;
TBparams.TBCatIVageBrac = [0,45; 45, Inf];
TBparams.prbCatIVrelapse =0.197;
%@using the AVERAGE parameters
%male   %female
TBparams.CatIVdeath = [...
0.016887412, 0.016887412;...
0.016887412, 0.016887412,...
];
%given not dead
TBparams.CatIVdefault = [...
0.01748523  ,   0.01748523  ;...
0.01748523  ,   0.01748523  ;...
];
%given have an outcome (ie, time to have outcome, not dead or defaulted)
TBparams.CatIVsuccess = [...
0.7375  ,   0.7375  ;...
0.7375  ,   0.7375  ;...
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sensitivity analysis: if transmission is more difficult for MDR
TBparams.MdrRelFitness = 0.7;  %base case used to be 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TBparams.modTransMatAgeBrac = [...
0   ,   4   ;...
5   ,   9   ;...
10  ,   14  ;...
15  ,   19  ;...
20  ,   24  ;...
25  ,   29  ;...
30  ,   34  ;...
35  ,   39  ;...
40  ,   44  ;...
45  ,   49  ;...
50  ,   54  ;...
55  ,   59  ;...
60  ,   64  ;...
65  ,   69  ;...
70  ,   Inf ;...
];
%Default: TBparams.contactScalingFactor = 1, riskInfectionPerContact = 6;
%10 infections per infected person a year
TBparams.riskInfectionPerContact = FOI_cal;
TBparams.contactScalingFactor = 1;
modContactMat = TBparams.contactScalingFactor*[...
0.69, 0.24, 0.36, 0.76, 0.79, 0.53, 0.97, 0.83, 0.38, 0.22, 0.36, 0.22, 0.22, 0.05, 0.09;...
0.59, 1.82, 0.43, 0.22, 0.79, 0.77, 0.42, 1.15, 1.09, 0.43, 0.21, 0.32, 0.15, 0.08, 0.09;...
0.25, 0.53, 1.69, 0.32, 0.09, 1.11, 0.69, 0.51, 0.58, 0.22, 0.14, 0.1, 0.2, 0.07, 0.1;...
0.18, 0.44, 0.79, 1.4, 0.38, 0.2, 0.67, 0.75, 0.47, 0.74, 0.47, 0.1, 0.1, 0.08, 0.08;...
0.42, 0.57, 0.17, 0.85, 1.02, 0.58, 0.31, 0.2, 0.27, 0.46, 0.22, 0.14, 0.1, 0.07, 0.19;...
0.61, 0.42, 0.44, 0.12, 0.49, 0.49, 0.42, 0.24, 0.29, 0.29, 0.17, 0.22, 0.14, 0.12, 0.1;...
0.57, 0.68, 0.32, 0.37, 0.28, 0.35, 0.8, 0.47, 0.25, 0.17, 0.15, 0.13, 0.15, 0.08, 0.05;...
0.74, 0.99, 0.51, 0.29, 0.46, 0.21, 0.21, 0.76, 0.49, 0.17, 0.14, 0.21, 0.2, 0.16, 0.13;...
0.18, 0.66, 0.69, 0.48, 0.23, 0.38, 0.27, 0.47, 0.45, 0.35, 0.19, 0.05, 0.19, 0.11, 0.24;...
0.2, 0.15, 0.27, 0.51, 0.44, 0.25, 0.38, 0.29, 0.4, 0.58, 0.27, 0.25, 0.13, 0.04, 0.29;...
0.36, 0.17, 0.2, 0.27, 0.2, 0.3, 0.23, 0.18, 0.12, 0.15, 0.33, 0.17, 0.09, 0.05, 0.23;...
0.15, 0.19, 0.19, 0.13, 0.28, 0.3, 0.2, 0.11, 0.19, 0.19, 0.35, 0.52, 0.24, 0.07, 0.22;...
0.18, 0.41, 0.17, 0.09, 0.11, 0.2, 0.23, 0.32, 0.23, 0.14, 0.2, 0.32, 0.27, 0.17, 0.11;...
0.26, 0.3, 0.26, 0.26, 0.15, 0.04, 0.15, 0.04, 0.15, 0.11, 0.04, 0.44, 0.3, 0.37, 0.22;...
0.07, 0.07, 0.17, 0.37, 0, 0.03, 0.07, 0.03, 0.3, 0.2, 0.07, 0.1, 0.13, 0.43, 0.7;...
    ]*TBparams.riskInfectionPerContact;
% modContactMat = TBparams.contactScalingFactor*[...
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0.8	0.5	0.5	0.5	0.5	0.8
% 0	0	0	0	0	0	0	0	0	0.5	0.8	0.5	0.5	0.5	0.5
% 0	0	0	0	0	0	0	0	0	0.5	0.5	0.8	0.5	0.5	0.5
% 0	0	0	0	0	0	0	0	0	0.5	0.5	0.5	0.8	0.5	0.5
% 0	0	0	0	0	0	0	0	0	0.8	0.5	0.5	0.5	0.8	0.5
% 0	0	0	0	0	0	0	0	0	0.5	0.8	0.5	0.5	0.5	0.8
%     ]*TBparams.riskInfectionPerContact;
TBparams.modMonRateMat = 30.5*modContactMat;

%only the the fraction of 0-15 year olds who are SS+ can transmit.
TBparams.fracKidsSSpos = 0.2;  %1 in the plos one version. articles 142 and 143.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TBparams.initialMdrProbInTBppl = 0.014796584;
TBparams.prbMDRtrtmtB4 = 0.197;
TBparams.prbMDRafterDefault = 0.242;
TBparams.prbMDRafterFailure = 0.187;
%seed the population with MDR in 1996
TBparams.MDRseedIn1996_lat = 0.005;   %covert latent.  number from the WHO_mdr_fracOfTotalPop.xls in the articles WHO fold.
TBparams.MDRseedIn1996_act = 0.005;   %covert active. number from the WHO_mdr_fracOfTotalPop.xls in the articles WHO fold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LATENT TO ACTIVE

TBparams.fastSlowActivationThreshold = 36;  %duration since infection (in months) before using slow activation probs.

%this is out_x from visualization.m in Dropbox\TBproject\Additional TB validation measure by age Feb 2014\ActivationArticles
TBparams.latentToActAgeBrac = [
0	4
5	9
10	14
15	19
20	24
25	29
30	34
35	39
40	44
45	49
50	54
55	59
60	64
65	69
70	74
75	79
80	84
85	89
90	94
95	Inf 
];

%TBsimParams for fast activation
%this out out_y from visualization.m
TBparams.timeSinceInf =  12*[0    0.5000    1.0000    1.5000    2.0000    2.5000    3.0000];  %months since infection

%vq from visualization.m
%  rows are ages and col are time from infection: 0, 0.5, 1, etc.
annualProbs = [
0.005638139	0.005638139	0.004936139	0.004197539	0.004197539	0.004116	0.004116
0.006891988	0.006891988	0.006006988	0.005121988	0.004714295	0.004116	0.004116
0.009632334	0.009632334	0.008747334	0.007476764	0.006116062	0.004116	0.004116
0.019881632	0.019881632	0.018611062	0.013870362	0.0118703	0.004116	0.004116
0.029022775	0.029022775	0.024282075	0.019429061	0.011674761	0.004116	0.004116
0.028190525	0.028190525	0.023449825	0.017473675	0.009719375	0.004116	0.004116
0.027358275	0.027358275	0.022617575	0.015518289	0.007763989	0.004116	0.004116
0.026526025	0.026526025	0.020774502	0.013562902	0.00756845	0.004116	0.004116
0.025693775	0.025693775	0.018819116	0.011607516	0.00756845	0.004116	0.004116
0.023873025	0.023873025	0.016661425	0.010428463	0.006976013	0.004116	0.004116
0.021628632	0.021628632	0.014417032	0.009582123	0.006129673	0.004116	0.004116
0.019384239	0.019384239	0.012172639	0.008735784	0.005283334	0.004116	0.004116
0.017139846	0.017139846	0.011186495	0.007889445	0.0051987	0.004116	0.004116
0.014895454	0.014895454	0.010340155	0.007043105	0.0051987	0.004116	0.004116
0.01271908	0.01271908	0.00942203	0.006242732	0.005160032	0.004116	0.004116
0.011154884	0.011154884	0.007857834	0.005856054	0.004773354	0.004116	0.004116
0.009590688	0.009590688	0.006293638	0.005469375	0.004386675	0.004116	0.004116
0.008026491	0.008026491	0.005082696	0.005082696	0.004116	0.004116	0.004116
0.006462295	0.006462295	0.004696018	0.004696018	0.004116	0.004116	0.004116
0.004898098	0.004898098	0.004309339	0.004309339	0.004116	0.004116	0.004116
];
annualRate = -log(1-annualProbs);
monthlyProbs_fastAct = 1-exp(-annualRate/12);  %monthly probs.  Need to scale this by the calibration factor.
%monthlyProbs_fastAct = 0.0035*ones(size(monthlyProbs_fastAct));  %%DEBUGGING
TBparams.latentToActLessT2Yr = latToAct_cal*monthlyProbs_fastAct;

%TBsimParams for slow activation
TBparams.timeSinceAct_slowBins = [1:0.5:100];  %years from activtion, with year three = 1 since the equation indexes at 1.
annualProbs = max(0.0001,-0.001*log(TBparams.timeSinceAct_slowBins)+ 0.0035);  %log version
%annualProbs = max(0.001,-0.001*log(TBparams.timeSinceAct_slowBins)+ 0.0035);  %DEBUGGING HIGHER FLOOR ON SLOW

%annualProbs = max(0.0001, 0.051129*exp(TBparams.timeSinceAct_slowBins+2)+ 0.00157);  %power version (plus 2 since third year from activation is TBparams.timeSinceAct_slowBins = 1

%annualProbs = 0.0035*ones(size(TBparams.timeSinceAct_slowBins)); %max(0.0001,-0.001*log(TBparams.timeSinceAct_slowBins)+ 0.0035);   %%DEBUGGING
monthlyProbs_slowAct = 1-exp(-1*( -log(1-annualProbs) ) /12);
TBparams.latentToActGreatT2Yr = latToAct_cal*monthlyProbs_slowAct;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DOTS RAMP UP
empiricalYear = 1996:1:2006;
%from wikipedia, double check once rntcp reports back up
empiricalCoverage = [...
0
0.075
0.15
0.225
0.3
0.4
0.5
0.730928222
0.923404649
0.97
1
];
coverageNeededAt = linspace(1996,2006,(12*10));
TBparams.DotsCoverageSequence = interp1(empiricalYear, empiricalCoverage, coverageNeededAt)';
TBparams.DotsFullCoverage = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of months needed for DST
TBparams.DSTperiodLength = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FULL DOTS PLUS RAMP UP SEQUENCE
TBparams.DotsPlusFullCoverage = 1;
%ramp up over 8 years (this is a 96 x 1 vector)  This was fitted to an exp
%using the rntcpt report dots plus data
TBparams.DotsPlusRampUpPeriod = 8*12;  %8 year ramp up
TBparams.DotsPlusRampUpSequence = 0.278.*exp(0.03109.*[1:1:96]'- 1.76.*ones(96,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TBparams.discountRate = 0.03;  %is this used anywhere??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%costs (from the valuesNeeded.xlsx in the CEA folder all in 2012USD)
%
TBparams.numTurnedAwayPerDiagnosed = 4.067841678;
TBparams.DSTspec = 1;
TBparams.DSTsensit = 1;

%RNTCP treatment
TBparams.monthlyNonmedPatientCosts = 6.327;
TBparams.SStestCost = 4.9264;
TBparams.DSTcost = 26.57162211;
TBparams.catIdrugCost = 3.599603981;
TBparams.catIclinicCost = 31.71427016;
TBparams.catIIdrugCost = 8.171993501;
TBparams.catIIclinicCost = 31.71427016;
TBparams.catIVdrugIPcost = 86.66;
TBparams.catIVdrugCPcost = 68.88;
TBparams.catIVclinicCost = 31.71427016;

%geneXpert
TBparams.geneXpertCost = 18.67629239;
TBparams.geneXpertVolDiscountCost = 14.94103391;
TBparams.geneXsensDS = 0.9040;
TBparams.geneXspecDS = 0.9840;
TBparams.geneXsensMDR = 0.9410;
TBparams.geneXspecMDR = 0.9700;

%private treatment
TBparams.monthlyPrivClinicCost = 69.55887983;
TBparams.PrivateTestCost = TBparams.SStestCost;  %default use SS for private diagnosis

%PPM
TBparams.perPatientPPMcost = 34.72322628;

%cost of living
TBparams.livingCostAgeBrac = [
0	9
10	19
20	29
30	39
40	49
50	59
60	69
70	Inf
];

TBparams.livingCost = [
0.204264181	0.187110683
0.171109992	0.261950115
1.229010851	2.633438033
1.618300179	2.779145581
1.655946108	2.92499098
2.527125792	2.366398641
3.313135235	2.440445391
4.756929941	4.178365436   
];
%smooth
medTime = mean(TBparams.livingCostAgeBrac,2);
medTime(find(medTime == Inf)) = 80;
TBparams.livingCost = interp1(medTime, TBparams.livingCost, [0:1:99]');
TBparams.livingCost(1:5,:) = repmat(TBparams.livingCost(6,:),5,1); 
TBparams.livingCost(81:100,:) = repmat(TBparams.livingCost(80,:),20,1); 

TBparams.livingCostAgeBrac = [[0:1:99]',[0:1:99]'];
TBparams.livingCostAgeBrac(end,end) = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%health state, treatment state, past treatment history must be <= value
TBparams.QalyHealthTreatState = [...
0   ,   0   ,  0 ;...
0   ,   0   ,  999;...
1   ,   0   ,  0 ;...
1   ,   0   ,  999 ;...
2   ,   0   ,  0 ;...
2   ,   0   ,  999 ;...
3   ,   0   ,  999 ;...
3   ,   1   ,  999 ;...
3   ,   2   ,  999 ;...
3   ,   3   ,  999 ;...
3   ,   4   ,  999 ;...
4   ,   0   ,  999 ;...
4   ,   1   ,  999 ;...
4   ,   2   ,  999 ;...
4   ,   3   ,  999 ;...
4   ,   4   ,  999 ;...
5   ,   0   ,  999 ;...
6   ,   0   ,  999 ;...
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for getting the life expectancies

TBparams.lebuilderCohortDurinYrs = 10;

i = 1;
for sex = 1:2
    for health = 0:4
        for trtStat = [0,1,2,4]
            for age = 0:TBparams.lebuilderCohortDurinYrs:99  
                if (health == 3 || health == 4)
                    PactTime = 1920;  %this is the startMonth jan 2026
                    PtimePerInfect = 1897;  %Feb 2024 (2 yrs from end simulation time)
                elseif (health == 1 || health == 2)
                    PactTime = 0;
                    PtimePerInfect = 1897;  %acquired disease 23 months ago
                else
                    PactTime = 0;
                    PtimePerInfect = 0;
                end
                LEbuildstate = [sex,age,0,0,2,0,health,trtStat,0,PtimePerInfect,0,0,0,PactTime,0,0,0,0,0,0];
                TBparams.LEbuilderCohort{i} = LEbuildstate;
                
                %out(i, 1:20) = LEbuildstate ;
                %out2(i, 1:4) = [sex, health, trtStat, age] ;
                i = i+1;
            end
        end
    end
end

%%copy and paste the resulting data into the excel in the outputs folder
% out = [[1:1:size(out,1)]',out];
% titleStr = ',Psex,Page,PageIncr,Psmoking,PurbanRur,PtrtmtCounter,Phealth,Ptreatment,PpastTreatment,PtimePerInfect,PtimeSmokingChange,PmdrTesting,PmdrTestingStatus,PactTime,PmdrEvolved,PseekTrtmt,PassignedSeekTrtmt,PtrtmtQ,PtrtmtQ_private,PprivTrtCount';
% tablePrinter(titleStr, out, 'LEbuilderStates', 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\');  %uncomment this, comment out the TBparams.LEbuilderCohort{i} = line, and remake if you change the cohorts
% tablePrinter('sex, health,trtStat, age', out2, 'LEbuilderStates_values', 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\');  %uncomment this, comment out the TBparams.LEbuilderCohort{i} = line, and remake if you change the cohorts


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%marks everyone for DST at the very beginning, for patients coming into
%normal treatment
TBparams.DSTinsteadOfSS = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for the minor revision PLOS ONE scenarios
TBparams.catIcuresMDR = 0;  %default 0
TBparams.catIIcuresMDR = 0.022;  %default 0
TBparams.MDRinDOTS2lessTransmissible = 0.022;  %default 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TBparams.limitOutputs = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Done running TBsimParams');
% 
% disp('NOW RUNNING SMOKING UNSTRATIFICATION');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Base mortality
% TBparams.fixedDeathYr = 0;
% 
% TBparams.baseMortAgeBrac = [...
% 0   ,   1   ;...
% 1   ,   4   ;...
% 5   ,   9   ;...
% 10  ,   14  ;...
% 15  ,   19  ;...
% 20  ,   24  ;...
% 25  ,   29  ;...
% 30  ,   34  ;...
% 35  ,   39  ;...
% 40  ,   44  ;...
% 45  ,   49  ;...
% 50  ,   54  ;...
% 55  ,   59  ;...
% 60  ,   64  ;...
% 65  ,   69  ;...
% 70  ,   74  ;...
% 75  ,   79  ;...
% 80  ,   84  ;...
% 85  ,   89  ;...
% 90  ,   94  ;...
% 95  ,   98  ;...
% 99  ,   Inf ];
% 
% datayears = [1990; 2000; 2009];
% 
% mortMaleSmokingOverTime = [...
% 0.007313128 0.005831265 0.004287449
% 0.000645625 0.000454064 0.000259966
% 0.000233306 0.000155821 0.00011416
% 0.000137491 0.000117493 8.08301E-05
% 0.000166653 0.000154155 0.000131658
% 0.000230807 0.000213311 0.00018415
% 0.000255801 0.000269964 0.00021831
% 0.00030412  0.000323281 0.000278295
% 0.000385759 0.000427409 0.000379095
% 0.000531525 0.000538188 0.000499042
% 0.000865459 0.000777198 0.000690595
% 0.001269194 0.001147674 0.00100283
% 0.001985526 0.001827495 0.001413999
% 0.003024583 0.00255257  0.002412086
% 0.004503993 0.004122313 0.003634214
% 0.006916801 0.005508937 0.005701186
% 0.0096333   0.008489591 0.008081336
% 0.013190565 0.010384868 0.011574165
% 0.018051746 0.013383796 0.01643921
% 0.024686837 0.018170391 0.023149508
% 0.033728086 0.025971778 0.032309847
% 1   1   1
% ];
% 
% mortFemaleSmokingOverTime=[...
% 0.007477735 0.005960499 0.004382037
% 0.000965367 0.000697257 0.000425743
% 0.000290791 0.00019748  0.00011916
% 0.000148322 0.000121659 8.49964E-05
% 0.000242471 0.000200813 0.000131658
% 0.000309952 0.000259133 0.000176651
% 0.000290791 0.00027163  0.000164986
% 0.000274962 0.000259133 0.000186649
% 0.000329112 0.00028246  0.000219976
% 0.000401586 0.000344107 0.000280794
% 0.000572336 0.000486548 0.000382427
% 0.000870454 0.000770536 0.000561509
% 0.001352418 0.001273355 0.00096953
% 0.002248302 0.00183914  0.001750965
% 0.003635874 0.003183256 0.002691372
% 0.005884286 0.004633399 0.004661601
% 0.007918482 0.007120363 0.006536874
% 0.012395864 0.009100834 0.010031842
% 0.018586762 0.012170335 0.015037458
% 0.026689859 0.017015244 0.022012444
% 0.036704442 0.024864003 0.031457098
% 1   1   1
% ];
% 
% mortMaleNonsmokingOverTime=[...
% 0.007313128 0.005831265 0.004287449
% 0.000645625 0.000454064 0.000259966
% 0.000233306 0.000155821 0.00011416
% 0.000137491 0.000117493 8.08301E-05
% 0.000166653 0.000154155 0.000131658
% 0.000230807 0.000213311 0.00018415
% 0.000255801 0.000269964 0.00021831
% 0.00030412  0.000323281 0.000278295
% 0.000385759 0.000427409 0.000379095
% 0.000531525 0.000538188 0.000499042
% 0.000865459 0.000777198 0.000690595
% 0.001269194 0.001147674 0.00100283
% 0.001985526 0.001827495 0.001413999
% 0.003024583 0.00255257  0.002412086
% 0.004503993 0.004122313 0.003634214
% 0.006916801 0.005508937 0.005701186
% 0.0096333   0.008489591 0.008081336
% 0.013190565 0.010384868 0.011574165
% 0.018051746 0.013383796 0.01643921
% 0.024686837 0.018170391 0.023149508
% 0.033728086 0.025971778 0.032309847
% 1   1   1
% ];
% 
% mortFemaleNonsmokingOverTime=[...
% 0.007477735 0.005960499 0.004382037
% 0.000965367 0.000697257 0.000425743
% 0.000290791 0.00019748  0.00011916
% 0.000148322 0.000121659 8.49964E-05
% 0.000242471 0.000200813 0.000131658
% 0.000309952 0.000259133 0.000176651
% 0.000290791 0.00027163  0.000164986
% 0.000274962 0.000259133 0.000186649
% 0.000329112 0.00028246  0.000219976
% 0.000401586 0.000344107 0.000280794
% 0.000572336 0.000486548 0.000382427
% 0.000870454 0.000770536 0.000561509
% 0.001352418 0.001273355 0.00096953
% 0.002248302 0.00183914  0.001750965
% 0.003635874 0.003183256 0.002691372
% 0.005884286 0.004633399 0.004661601
% 0.007918482 0.007120363 0.006536874
% 0.012395864 0.009100834 0.010031842
% 0.018586762 0.012170335 0.015037458
% 0.026689859 0.017015244 0.022012444
% 0.036704442 0.024864003 0.031457098
% 1   1   1
% ];
% 
% 
% mortMaleSmoking = interp1(datayears, mortMaleSmokingOverTime', [1990:1:2009]);
% mortFemaleSmoking = interp1(datayears, mortFemaleSmokingOverTime', [1990:1:2009]);
% mortMaleNonsmoking = interp1(datayears, mortMaleNonsmokingOverTime', [1990:1:2009]);
% mortFemaleNonsmoking = interp1(datayears, mortFemaleNonsmokingOverTime', [1990:1:2009]);
% 
% for i = 1 : ( max(datayears) - min(datayears) + 1)
%     TBparams.mortSmoking{i} = [mortMaleSmoking(i,:)', mortFemaleSmoking(i,:)'];
%     TBparams.mortNonsmoking{i} = [mortMaleNonsmoking(i,:)', mortFemaleNonsmoking(i,:)'];
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% smokingPrev = [... 
% 0.0000000   0.0000000    
% 0.0000000   0.0000000
% 0.0000000   0.0000000
% 0.0000000   0.0000000
% 0.0540131   0.0025893
% 0.1617864   0.0068130
% 0.2213199   0.0141294
% 0.2808533   0.0214457
% 0.3190453   0.0301124
% 0.3572372   0.0387790
% 0.3712019   0.0489459
% 0.3851665   0.0591128
% 0.3597671   0.0727734
% 0.3343677   0.0864339
% 0.3258945   0.0958011
% 0.3174212   0.1051682
% 0.2916656   0.0959240
% 0.2574367   0.0960470
% 0.2627954   0.1343809
% 0.2681540   0.1727147
% 0.1457413   0.1142369
% 0.0233286   0.0557590
% ];  %these are the GATS numbers that Jeremy sent me.   
% %fixed so that age brackets are the same as mortality tables
% 
% maleFemaleLeftAlive = [... 
% 100000  100000
% 93311   93169
% 91303   90112
% 90455   89049
% 89820   88401
% 88993   87342
% 87862   85993
% 86449   84601
% 84787   83297
% 82639   81899
% 80012   80225
% 76366   77917
% 71280   74395
% 63865   68918
% 54768   61703
% 42691   50930
% 30552   38475
% 18106   24891
% 9470    14176
% 4018    6559
% 1355    2384
% 339 634
% ];  %at every age.  Info from 2000 WHO lifetables
% 
% unsexedSmokingPrev = sum((smokingPrev.*maleFemaleLeftAlive),2) ./ sum(maleFemaleLeftAlive,2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %mortality prob from non-treated TB;
% untreatedTBmortProbNonSmoker = 0.02469;  
% smokingTBmortRiskRatio = 1.475;
% 
% untreatedTBmortProbSmoker = untreatedTBmortProbNonSmoker * smokingTBmortRiskRatio;
% 
% warning off
% untreatedTBmortRateNonSmok = -log(1-untreatedTBmortProbNonSmoker);
% untreatedTBmortRateSmok = -log(1-untreatedTBmortProbSmoker);
% warning on
% 
% 
% for i = 1: ( max(datayears) - min(datayears) + 1)
%     warning off
%     baseMortSmokingRate{i} = -log(1-TBparams.mortSmoking{i});
%     baseMortNonsmokingRate{i} = -log(1-TBparams.mortNonsmoking{i});
%     warning on
% 
%     TBparams.nonDOTsTBmortSmoking{i} = 1-exp(-(baseMortSmokingRate{i} + untreatedTBmortRateSmok));
%     TBparams.nonDOTsTBmortNonSmok{i} = 1-exp(-(baseMortNonsmokingRate{i} + untreatedTBmortRateNonSmok));
%     
%     remixedNonDOTsTBmortSmoking{i} = (smokingPrev).*TBparams.nonDOTsTBmortSmoking{i} + (1-smokingPrev).*TBparams.nonDOTsTBmortNonSmok{i};
%     TBparams.nonDOTsTBmortSmoking{i} = remixedNonDOTsTBmortSmoking{i};
%     TBparams.nonDOTsTBmortNonSmok{i} = remixedNonDOTsTBmortSmoking{i};
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %LATENT TO ACTIVE
% 
% %TBparams.latentToActiveTB = 0.000941662;
% TBparams.latentToActAgeBrac = [...
% 0   ,   1   ;...
% 1   ,   4   ;...
% 5   ,   9   ;...
% 10  ,   14  ;...
% 15  ,   19  ;...
% 20  ,   24  ;...
% 25  ,   29  ;...
% 30  ,   34  ;...
% 35  ,   39  ;...
% 40  ,   44  ;...
% 45  ,   49  ;...
% 50  ,   54  ;...
% 55  ,   59  ;...
% 60  ,   64  ;...
% 65  ,   69  ;...
% 70  ,   74  ;...
% 75  ,   79  ;...
% 80  ,   84  ;...
% 85  ,   89  ;...
% 90  ,   94  ;...
% 95  ,   98  ;...
% 99  ,   Inf ];
% 
% TBparams.latentToActLessT2Yr = [...
% 0.000449899
% 0.000449899
% 0.000449899
% 0.000099995
% 0.000099995
% 0.000466558
% 0.000466558
% 0.000466558
% 0.000349939
% 0.000349939
% 0.000349939
% 0.000349939
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% 0.000141657
% ];
% 
% TBparams.latentToActGreatT2Yr = [...
% 0.00019998
% 0.00019998
% 0.00019998
% 0.00011666
% 1.00011666
% 0.000158321
% 0.000158321
% 0.000158321
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% 0.000099995
% ];
% 
% smokerActivationRiskRatio = 2.3;
% 
% TBparams.latentToActLessT2YrSmoker = TBparams.latentToActLessT2Yr * smokerActivationRiskRatio;
% TBparams.latentToActGreatT2YrSmoker = TBparams.latentToActGreatT2Yr * smokerActivationRiskRatio;
% 
% remixedActivationLess2Yr = (unsexedSmokingPrev).*TBparams.latentToActLessT2YrSmoker + (1-unsexedSmokingPrev).*TBparams.latentToActLessT2Yr;
% remixedActivationGreat2Yr = (unsexedSmokingPrev).*TBparams.latentToActGreatT2YrSmoker + (1-unsexedSmokingPrev).*TBparams.latentToActGreatT2Yr;
% 
% TBparams.latentToActLessT2Yr = remixedActivationLess2Yr;
% TBparams.latentToActLessT2YrSmoker = remixedActivationLess2Yr;
% TBparams.latentToActGreatT2Yr = remixedActivationGreat2Yr;
% TBparams.latentToActGreatT2YrSmoker = remixedActivationGreat2Yr;
% 
% disp('DONE RUNNING SMOKING UNSTRATIFICATION');

