
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Base mortality
TBparams.fixedDeathYr = 0;

TBparams.baseMortAgeBrac = [...
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
95  ,   98  ;...
99  ,   Inf ];

datayears = [1990; 2000; 2009];

mortMaleSmokingOverTime = [...
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
1	1	1
];

mortFemaleSmokingOverTime=[...
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
1	1	1
];

mortMaleNonsmokingOverTime=[...
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
1	1	1
];

mortFemaleNonsmokingOverTime=[...
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
1	1	1
];


mortMaleSmoking = interp1(datayears, mortMaleSmokingOverTime', [1990:1:2009]);
mortFemaleSmoking = interp1(datayears, mortFemaleSmokingOverTime', [1990:1:2009]);
mortMaleNonsmoking = interp1(datayears, mortMaleNonsmokingOverTime', [1990:1:2009]);
mortFemaleNonsmoking = interp1(datayears, mortFemaleNonsmokingOverTime', [1990:1:2009]);

for i = 1 : ( max(datayears) - min(datayears) + 1)
	TBparams.mortSmoking{i} = [mortMaleSmoking(i,:)', mortFemaleSmoking(i,:)'];
	TBparams.mortNonsmoking{i} = [mortMaleNonsmoking(i,:)', mortFemaleNonsmoking(i,:)'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smokingPrev = [... 
0.0000000	0.0000000    
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
];  %these are the GATS numbers that Jeremy sent me.   
%fixed so that age brackets are the same as mortality tables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mortality prob from non-treated TB;
untreatedTBmortProbNonSmoker = 0.02469;  
smokingTBmortRiskRatio = 1.475;

untreatedTBmortProbSmoker = untreatedTBmortProbNonSmoker * smokingTBmortRiskRatio;

warning off
untreatedTBmortRateNonSmok = -log(1-untreatedTBmortProbNonSmoker);
untreatedTBmortRateSmok = -log(1-untreatedTBmortProbSmoker);
warning on


for i = 1: ( max(datayears) - min(datayears) + 1)
    warning off
    baseMortSmokingRate{i} = -log(1-TBparams.mortSmoking{i});
    baseMortNonsmokingRate{i} = -log(1-TBparams.mortNonsmoking{i});
    warning on

    TBparams.nonDOTsTBmortSmoking{i} = 1-exp(-(baseMortSmokingRate{i} + untreatedTBmortRateSmok));
    TBparams.nonDOTsTBmortNonSmok{i} = 1-exp(-(baseMortNonsmokingRate{i} + untreatedTBmortRateNonSmok));
    
    remixedNonDOTsTBmortSmoking{i} = (smokingPrev).*TBparams.nonDOTsTBmortSmoking{i} + (1-smokingPrev).*TBparams.nonDOTsTBmortNonSmok{i};
    TBparams.nonDOTsTBmortSmoking{i} = remixedNonDOTsTBmortSmoking{i};
    TBparams.nonDOTsTBmortNonSmok{i} = remixedNonDOTsTBmortSmoking{i};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LATENT TO ACTIVE

%TBparams.latentToActiveTB = 0.000941662;
TBparams.latentToActAgeBrac = [...
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
95  ,   98  ;...
99  ,   Inf ];

TBparams.latentToActLessT2Yr = [...
0.000449899
0.000449899
0.000449899
0.000099995
0.000099995
0.000466558
0.000466558
0.000466558
0.000349939
0.000349939
0.000349939
0.000349939
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
0.000141657
];

TBparams.latentToActGreatT2Yr = [...
0.00019998
0.00019998
0.00019998
0.00011666
1.00011666
0.000158321
0.000158321
0.000158321
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
0.000099995
];

smokerActivationRiskRatio = 2.3;

TBparams.latentToActLessT2YrSmoker = TBparams.latentToActLessT2Yr * smokerActivationRiskRatio;
TBparams.latentToActGreatT2YrSmoker = TBparams.latentToActGreatT2Yr * smokerActivationRiskRatio;

remixedActivationLess2Yr = (smokingPrev).*TBparams.latentToActLessT2YrSmoker + (1-smokingPrev).*TBparams.latentToActLessT2Yr;
remixedActivationGreat2Yr = (smokingPrev).*TBparams.latentToActGreatT2YrSmoker + (1-smokingPrev).*TBparams.latentToActGreatT2Yr;

TBparams.latentToActLessT2Yr = remixedActivationLess2Yr;
TBparams.latentToActLessT2YrSmoker = remixedActivationLess2Yr;
TBparams.latentToActGreatT2Yr = remixedActivationGreat2Yr;
TBparams.latentToActGreatT2YrSmoker = remixedActivationGreat2Yr;
