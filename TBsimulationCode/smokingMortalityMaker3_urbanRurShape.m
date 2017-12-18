function output = smokingMortalityMaker3_urbanRurShape

%calculation of smoking mortality
%matracies calculated on excel spreadsheet

t = 1/12;

ageBrac = [...
20	,	24	;...
25	,	29	;...
30	,	34	;...
35	,	39	;...
40	,	44	;...
45	,	49	;...
50	,	54	;...
55	,	59	;...
60	,	64	;...
65	,	69	;...
70	,	74	;...
];

%for all the following matracies, the first column is men, the 2nd women.
%men        %women

%from Indian smoking prevalence data from Jeremy
fracSmoker = [...
0.144027	,	0.020664	;...
0.302404	,	0.027755	;...
0.381743	,	0.038614	;...
0.415032	,	0.049292	;...
0.430601	,	0.036372	;...
0.424791	,	0.05139	    ;...
0.413023	,	0.070911	;...
0.406085	,	0.054887	;...
0.381571	,	0.059954	;...
0.300339	,	0.047119	;...
0.284613	,	0.0416405	;...
];

%this is the difference in death rate for smokers and nonsmokers(from
%paper)
riskRatio = [...
1.2	,	1.5	;...
1.2	,	1.5	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.7	,	2	;...
1.6	,	1.3	;...
];

%these are death rates for the whole Indian population, taken from WHO
%lifetables
PrOverallDeath = [...
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
];

fracNonSmokers = 1-fracSmoker;
PrNonSmokerDeath = PrOverallDeath ./ ((fracSmoker.*riskRatio) + fracNonSmokers);
PrSmokerDeath = riskRatio .* PrNonSmokerDeath;

PrSmokerDeath = [...
0.000248806	,	0.000384725	;...
0.000305481	,	0.000401868	;...
0.000385609	,	0.000357987	;...
0.000572696	,	0.00049651	;...
0.000714668	,	0.000611168	;...
0.00103288	,	0.000849286	;...
0.001535571	,	0.001324673	;...
0.002369741	,	0.002068539	;...
0.003350996	,	0.002977676	;...
0.00544979	,	0.004923314	;...
0.007528648	,	0.005949102	;...
];

PrNonSmokerDeath = [...
0.000207338	,	0.000256483	;...
0.000254567	,	0.000267912	;...
0.000284797	,	0.000255163	;...
0.000324328	,	0.000271362	;...
0.000404728	,	0.000334027	;...
0.000588377	,	0.000466897	;...
0.000874733	,	0.000728243	;...
0.001456738	,	0.001227175	;...
0.00205994	,	0.001766527	;...
0.003552476	,	0.003097212	;...
0.004705405	,	0.004576232	;...
];


%folderName = 'C:\MATLAB701\work\smokingMortality';
folderName = 'C:\Users\Sze\Dropbox\TBproject\code\smokingMortality';

output = {PrNonSmokerDeath, PrSmokerDeath};


%other Jha paper for errors (article 2)
otherJhaAgeBrac = [...
30;...
40;...
50;...
60;...
67.5;...
72.5;...
75;...
];

urbanAve = [...
2.4	;...
2.7	;...
2.4	;...
1.9	;...
1.8	;...
1.4	;...
1.3	;...
];

urbanLower =[...
2.1	;...
2.4	;...
2.2	;...
1.7	;...
1.6	;...
1.2	;...
1.1	;...
];

urbanUpper =[...
2.9	;...
3	;...
2.6	;...
2.1	;...
2	;...
1.6	;...
1.4	;...
];

ruralAve = [...
1.1	;...
1.6	;...
1.7	;...
1.7	;...
1.6	;...
1.2	;...
1.3	;...
];

ruralLower =[...
1	;...
1.4	;...
1.6	;...
1.5	;...
1.3	;...
1	;...
1.1	;...
];

ruralUpper =[...
1.4	;...
1.8	;...
2	;...
1.9	;...
1.9	;...
1.5	;...
1.6	;...
];

riskRatio = output{2}./output{1};
plot(ageBrac(:,1),riskRatio);
title('Risk Ratios over age');
xlabel('age');
legend('men','women');
print('-dpng','-r300',[folderName '/riskRatios.png']);
hold on
errorbar(otherJhaAgeBrac, urbanAve, urbanAve-urbanLower, urbanUpper-urbanAve,'ko');
hold on
errorbar(otherJhaAgeBrac, ruralAve, ruralAve-ruralLower, ruralUpper-ruralAve,'ro');
title('Risk Ratios over age, errors from other Jha paper (black urban, red rural)');
print('-dpng','-r300',[folderName '/riskRatios_withErrors.png']);


figure
outputAll = [output{1},output{2}];
plot(ageBrac(:,1),outputAll);
leg = legend('nonSmoking men', 'nonSmoking women', 'smoking men', 'smoking women');
set(leg,'Location','NorthWest');
title('Probabilities of Death');
xlabel('age');
print('-dpng','-r300',[folderName '/probDeath_constantRR.png']);

figure
diff = output{2} - output{1};
plot(ageBrac(:,1),diff);
leg = legend('men','women');
set(leg,'Location','NorthWest');
title('difference in death prob');
xlabel('age');
print('-dpng','-r300',[folderName '/diff_constantRR.png']);

%this is the difference in death rate for smokers and nonsmokers(from
%paper)
USpaperAgeBrac = [...
30	,	34	;...
35	,	39	;...
40	,	44	;...
45	,	49	;...
50	,	54	;...
55	,	59	;...
60	,	64	;...
65	,	69	;...
70	,	74	;...
];

USpaperDiff = [...
2.06404;...
2.42129;...
2.38723;...
2.13251;...
1.81239;...
1.51751;...
1.28242;...
1.14802;...
1.01988;...
];
figure
plot(USpaperAgeBrac(:,1),USpaperDiff);
title('risk ratio from US paper (men only)');
xlabel('age');
print('-dpng','-r300',[folderName '/US_riskRatios.png']);
%



