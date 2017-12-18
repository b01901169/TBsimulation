%TBparams Overwriter
%each cluster should only have one thing turned on (set to 1).
%forceOfInfectionParam and latToActiveScalingFactor are not dummies.

%DOTS treatment quality
noTreat = 0; %default 0
worstTreat = 0; %default 0  WRONG NOW IF YOU NEED THIS LOOK AT DETAILEDANALYSISPLAN_DATATABLE_SUMMER3 TREATMENT
aveTreat = 1; %default 1
bestTreat = 0; %default 0  WRONG NOW IF YOU NEED THIS LOOK AT DETAILEDANALYSISPLAN_DATATABLE_SUMMER3 TREATMENT
perfectTreat = 0; %default 0

%DOTS uptake
womUptake = 0; %default 0
empUptake = 1; %default 1.
ultraUptake = 0; %default 0.

%DOTS+ Ramp up schedule
noDotsPlus = 0; %default 0
halfDotsPlus = 0; %default 0
fullDotsPlus = 1;  %default 1

%DOTS+ quality
perfectDotsPlus = 0; %default 0

%background params
calibrationRun = 0;  %this is if you want the 230 with the life expectancy, etc.
calibrationRunYr = 1990;  %this can be 1990 or 2000.  only works if calibrationRun == 1

%also for calibrationRuns make sure to turn treatment off
noSmoke = 0;  %default 0

%Actual GeneXpert CEA Scenarios (all default 0)
GeneXatInitialRNTCP = 0;
GeneXatDstRNTCP = 0;

%PPM CEA Scneario
PPMlevel = 0;  %default 0, can range from 0 to 1.
TBparams.PPMeffectiveness = 1; %default 1, can range from 0 to 1.

%GeneXpert Scenarios (all default 0)
GeneXinsteadOfSS = 0;  %perfect diagnosis of DS TB (scenario c)
GeneXinsteadOfLJ = 0; %fast diagnosis on MDR TB (scenario b)
GeneXinitialDST = 0; %fast diagnosis of MDR TB at the first clinic visit (scenario d)


%%TB MAC scenarios
TBparams.TBmacMultiplier =1;
TBMAC1a_increaseTrt = 0; 
TBMAC1b_increaseTrt = 0;
TBMAC2a_initDefReduction = 0;  
TBMAC2b_trtSuccess = 0; 
TBMAC2c_MDRtrtSuccess = 0; 
TBMAC3_XpertReplacesSmear = 0;    
TBMAC3b_XpertReplacesSmear_SSsensUp = 0; 
TBMAC4_activeCaseFinding = 0;

% the analogous advocate scenarios
TBMAC1a_adv = 0;
TBMAC1b_adv = 0;
TBMAC2a_adv = 0; 
TBMAC2b_adv = 0;
TBMAC2c_adv = 0;
TBMAC3a_adv = 0; 
TBMAC3b_adv = 0;
Tbmac4_adv = 0;
TBMAC5_preventativeTherapyForLat = 0; %default 0, can be changed to 1

%other
otherSpecialParamRuns = 0;  %default 0.  This section can be changed to whatever the particular run calls for.  At the very bottom.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%SMOKING PREVALENCE%%%%%%%%%%%%%
%Default: smoking on
if noSmoke == 1
    %NO SMOKING AT ALL
    TBparams.smokingChurnAgeBrac = [5:1:100;6:1:101]';
    TBparams.smoker_quit_ratio = 1;
    TBparams.nonsmoker_start_ratio = 0;
    TBparams.smokingPrev = zeros(11,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%case detection probability for sputum smear process;
%Default: treatment on
if noTreat == 1
    %TURN TREATMENT OFF FOR DEBUGGING PURPOSES, DELETE THIS IF DONE DEBUGGING
    % %NO TREATMENT
    TBparams.SSsensit = 0; %SS test has 50% to 70% sensitivity; base is 60%
    TBparams.rapidDiagSens = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default: Average
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
if worstTreat == 1
    %@@@@@@@@@@@using the WORST parameters
    %male   %female
    TBparams.CatIdeath = repmat([ 0.04645828  ,   0.04645828 ],size(TBparams.TBCatIageBrac,1),1);
    %given not dead
    TBparams.CatIdefault = [...
        0.107110335 0.108790199
        0.105012557 0.080001751
        0.096394837 0.079411698
        0.07475726  0.068831005
        0.102370787 0.055991854
        0.09984767  0.051735193
        0.10137445  0.091578793
        ];
    %given have an outcome (ie, time to have outcome, not dead or defaulted)
    TBparams.CatIsuccess = repmat([0.90 ,   0.90],size(TBparams.TBCatIageBrac,1),1);
    
    %CAT II
    TBparams.CatIIdeath = repmat([ 0.04135318  ,   0.04135318 ],size(TBparams.TBCatIageBrac,1),1);
    TBparams.CatIIdefault = [...
        0.26637387  0.270551542
        0.261156882 0.198957233
        0.23972538  0.197489826
        0.185914651 0.171176584
        0.254587036 0.139246756
        0.248312267 0.128660819
        0.252109233 0.227748307
        ];
    TBparams.CatIIsuccess = repmat([0.83    ,   0.83],size(TBparams.TBCatIageBrac,1),1);
    
end
if bestTreat == 1
    %@@@@@@@@@@@using the BEST parameters
    TBparams.CatIdeath = repmat([ 0.00267381    ,   0.00267381 ],size(TBparams.TBCatIageBrac,1),1);
    %given not dead
    TBparams.CatIdefault = [...
        0.005103163 0.005183198
        0.005003216 0.003811602
        0.004592634 0.003783489
        0.003561733 0.003279383
        0.004877352 0.002667675
        0.00475714  0.00246487
        0.004829882 0.004363178
        ];
    %given have an outcome (ie, time to have outcome, not dead or defaulted)
    TBparams.CatIsuccess = repmat([0.99 ,   0.99],size(TBparams.TBCatIageBrac,1),1);
    
    %CAT II
    TBparams.CatIIdeath = repmat([ 0.014543831  ,   0.014543831],size(TBparams.TBCatIageBrac,1),1);
    TBparams.CatIIdefault = [...
0.00827161  0.008401337
0.008109608 0.006178146
0.007444104 0.006132579
0.005773139 0.005315483
0.007905597 0.004323978
0.007710749 0.003995257
0.007828655 0.007072184
        ];
    TBparams.CatIIsuccess = repmat([0.97    ,   0.97],size(TBparams.TBCatIageBrac,1),1);
end
if perfectTreat == 1
    %@@@@@@@@@@@using the PERFECT parameters
    %male   %female
    TBparams.CatIdeath = repmat([ 0 ,   0 ],size(TBparams.TBCatIageBrac,1),1);
    %given not dead
    TBparams.CatIdefault = repmat([ 0   ,   0 ],size(TBparams.TBCatIageBrac,1),1);    
    %given have an outcome (ie, time to have outcome, not dead or defaulted)
    TBparams.CatIsuccess = repmat([1,1],size(TBparams.TBCatIageBrac,1),1);   
    
    %CAT II
    TBparams.CatIIdeath = repmat([ 0    ,   0 ],size(TBparams.TBCatIageBrac,1),1);
    %given not dead
    TBparams.CatIIdefault = repmat([ 0  ,   0 ],size(TBparams.TBCatIageBrac,1),1);    
    %given have an outcome (ie, time to have outcome, not dead or defaulted)
    TBparams.CatIIsuccess = repmat([1,1],size(TBparams.TBCatIageBrac,1),1);   
end
if aveTreat == 1
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%CEA and REAL GeneXpert Scenarios%%%%%%%%%%%%%%%%%%%%%%%

if GeneXatInitialRNTCP == 1
%geneX at the beginning of RNTCP
TBparams.DSTinsteadOfSS  = 1;  %do dst after initial entry
TBparams.DSTperiodLength = 0;
TBparams.SSsensit = TBparams.geneXsensDS; 
TBparams.SSspec = TBparams.geneXspecDS;
TBparams.DSTspec = TBparams.geneXspecMDR;
TBparams.DSTsensit = TBparams.geneXsensMDR;
TBparams.DSTcost = TBparams.geneXpertCost ;  %can be TBparams.geneXpertVolDiscountCost
TBparams.SStestCost = TBparams.geneXpertCost ;
end

if GeneXatDstRNTCP  == 1
%geneX for DST at the current treatment algorithm only
TBparams.DSTperiodLength = 0;
TBparams.DSTspec = TBparams.geneXspecMDR;
TBparams.DSTsensit = TBparams.geneXsensMDR;
TBparams.DSTcost = TBparams.geneXpertCost ;  %can be TBparams.geneXpertVolDiscountCost
end

if PPMlevel > 0
%PPM
TBparams.privateReferToRNTCP = PPMlevel ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%OLD GENEXPERT SCENARIOS %%%%%%%%%%%%%%%%%
if GeneXinsteadOfSS == 1
    TBparams.SSsensit = 1;
end
if GeneXinsteadOfLJ == 1
    TBparams.DSTperiodLength = 0;
end
if GeneXinitialDST == 1
    TBparams.DSTinsteadOfSS  = 1;
    TBparams.DSTperiodLength = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%UPTAKE%%%%%
%knowledge of tb for treatment uptake
%WOMEN AS GOOD AS MEN
if womUptake == 1
    
end
if empUptake == 1
    
end
if ultraUptake == 1
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%perfectDotsPlus
if perfectDotsPlus == 1
    %Cat IV death prob;
    
    TBparams.TBCatIVageBrac = [0,45; 45, Inf];
    
    TBparams.prbCatIVrelapse =0.197;
    
    %@using the AVERAGE parameters
    %male   %female
    TBparams.CatIVdeath = [...
        0, 0;...
        0, 0,...
        ];
    
    %given not dead
    TBparams.CatIVdefault = [...
        0  ,   0  ;...
        0  ,  0  ;...
        ];
    
    %given have an outcome (ie, time to have outcome, not dead or defaulted)
    TBparams.CatIVsuccess = [...
        1  ,  1  ;...
        1  ,  1  ;...
        ];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DOTS PLUS RAMP UP SEQUENCE
if noDotsPlus == 1
    %none - no dots plus at all
    TBparams.DotsPlusFullCoverage = 0;
    TBparams.DotsPlusRampUpSequence = zeros(96, 1);
end
if halfDotsPlus == 1
    %HALF DOTS PLUS RAMP UP SEQUENCE
    TBparams.DotsPlusFullCoverage = TBparams.DotsPlusRampUpSequence(48,1);
    %ramp up over 8 years (this is a 96 x 1 vector)  This was fitted to an exp
    %using the rntcpt report dots plus data
    TBparams.DotsPlusRampUpSequence = 0.278.*exp(0.03109.*[1:1:96]' - 1.76.*ones(96,1));
    TBparams.DotsPlusRampUpSequence(49:end,1) = TBparams.DotsPlusRampUpSequence(48,1);
end
if fullDotsPlus == 1
    %FULL DOTS PLUS RAMP UP SEQUENCE
    TBparams.DotsPlusFullCoverage = 1;
    %ramp up over 8 years (this is a 96 x 1 vector)  This was fitted to an exp
    %using the rntcpt report dots plus data
    TBparams.DotsPlusRampUpSequence = 0.278.*exp(0.03109.*[1:1:96]' - 1.76.*ones(96,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calibration run
if calibrationRun == 1
    TBparams.MDRseedIn1996 = 0;
    
    TBparams.hisPopGrowthPerc = [...
        0.02076
        0.02076
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        0
        ];
    
    TBparams.constantHisPopGrowth = 0.014;
    
    if calibrationRunYr == 1990
        TBparams.fixedDeathYr = 1;
    elseif calibrationRunYr == 2000
        TBparams.fixedDeathYr = 11;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TB Mac scenario 1
if TBMAC1a_increaseTrt == 1;     %increase people accessing trt (jan 2016 to dec 2022)
    TBmac1ScalingFac = [1:(6/83):7]*TBparams.TBmacMultiplier;
    TBmac1ScalingFac = [1:(1/83):2]*TBparams.TBmacMultiplier;
end
if TBMAC1b_increaseTrt == 1  %increase proportion in high quality care (jan 2016 to dec 2022)    
    TBparams.PPMeffectivenessVal = [0:(0.8/83):0.8] *TBparams.TBmacMultiplier;
end

%TB Mac scenario 2 
if TBMAC2a_initDefReduction == 1   %2a (implemented after b and c so effects will be on top if both done)
    TbMac2_initDefaultDec = 0.5*ones(1,84)*TBparams.TBmacMultiplier;  %instantaneous
    TbMac2_MDRinitDefaultDec = [[1:-((6/11)/59):(5/11)], (5/11)*ones(1,83-59)]*TBparams.TBmacMultiplier;  %(jan 2016 to dec 2020)
end
if TBMAC2b_trtSuccess == 1  %%2b (jan 2016 to dec 2022)
    cat12initScalar = [1:-(0.5/83):0.5]*TBparams.TBmacMultiplier;  %scalar on initial default in cat 1 and 2
    cat12DeathInitScalar = [1:-((5/11)/83):(6/11)]*TBparams.TBmacMultiplier; 
end
if TBMAC2c_MDRtrtSuccess == 1   %%2c (jan 2016 to dec 2022)
    cat4initScalar = [1:-(0.5/83):0.5]*TBparams.TBmacMultiplier;  %scalar on default in cat4
    cat4deathScalar = [1:-((10/22)/83):(12/22)]*TBparams.TBmacMultiplier;  %scalar on cat4 death
end

%TB Mac scenario3
if (TBMAC3_XpertReplacesSmear == 1 || TBMAC3b_XpertReplacesSmear_SSsensUp == 1)
    GeneXcoverageMat = [0:(0.3/47):0.3]*TBparams.TBmacMultiplier;
end

%TB Mac scenario4 and 5
if (TBMAC4_activeCaseFinding == 1)
    TBmac4_sens = 0.90; TBmac4_spec = 0.9975; TBmac4_freq = 12; ACFpercOfUntrted = 0.016*TBparams.TBmacMultiplier;
end

%% advocate numbers -- make sure the vector lengths line up with above and you hit the trigger

if (TBMAC1a_adv == 1)
    TBMAC1a_increaseTrt = 1;
    TBmac1ScalingFac = [[1:(14/59):15], 15*ones(1,83-59)];
end
if TBMAC1b_adv == 1  %increase proportion in high quality care (jan 2016 to dec 2022) 
    TBMAC1b_increaseTrt = 1;
    TBparams.PPMeffectivenessVal = [[0:(0.95/59):0.95], 0.95*ones(1,83-59)];
end

if TBMAC2a_adv == 1 %2a_adv
    TBMAC2a_initDefReduction = 1;
    TbMac2_initDefaultDec = [[1:-(1/59):0],zeros(1,83-59)];  
    TbMac2_MDRinitDefaultDec = [[1:-(1/59):0],zeros(1,83-59)];  %(jan 2016 to dec 2020)
end
if TBMAC2b_adv == 1
    TBMAC2b_trtSuccess = 1;
    cat12initScalar = [[1:-(0.75/59):0.25],0.25*ones(1,83-59)];  %scalar on initial default in cat 1 and 2
    cat12DeathInitScalar = [[1:-(0.73/59):(0.27)],0.27*ones(1,83-59)]; 
end    
if  TBMAC2c_adv == 1   %2c_adv
    TBMAC2c_MDRtrtSuccess = 1;   %%2c (jan 2016 to dec 2022)
    cat4initScalar = [[1:-(0.75/59):0.25],0.25*ones(1,83-59)];  %scalar on default in cat4
    cat4deathScalar = [[1:-((16/22)/59):(6/22)],(6/22)*ones(1,83-59)];  %scalar on cat4 death
    cat4failureScalar = [[0:(1.2/59):1.2],1.2*ones(1,83-59)];   
end

if TBMAC3a_adv == 1
    TBMAC3_XpertReplacesSmear = 1;
    GeneXcoverageMat = [0:(1/47):1];
end
if TBMAC3b_adv == 1
    TBMAC3b_XpertReplacesSmear_SSsensUp = 1;
    GeneXcoverageMat = [0:(1/47):1];
end

if (Tbmac4_adv == 1 || TBMAC5_preventativeTherapyForLat == 1)
    TBMAC4_activeCaseFinding = 1;
    Tbmac4_adv = 1;
    TBmac4_sens = 0.90; TBmac4_spec = 0.9975; TBmac4_freq = 6; ACFpercOfUntrted_vec = [[0:0.3/11:0.3],0.3*ones(1,11)];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if otherSpecialParamRuns == 1
    %put whatever special things you want in here

end
