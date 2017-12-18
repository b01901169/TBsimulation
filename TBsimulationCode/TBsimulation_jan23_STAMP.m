function TBsimulation_jan23(folderName, logComment,durationYrs, numberPpl, plotResolution, loadBurnInStr, startScenarioYr,startScenarioYr2, latToAct_cal, oldActSlope_cal, oldActIntercept_cal, FOI_cal, aveUptake_cal, cat2uptake_cal, simParamsFolder)
% Name: TBsimulation.m
% Date: May 02, 2011, radically revised July 6, 2011
% Most Recently Updated: June 5, 2014
% Author: Sze Suen
% Function Call:
%function TBsimulation_july6(folderName, logComment,durationYrs, numberPpl, plotResolution, loadBurnInStr, startScenarioYr,startScenarioYr2, latToAct_cal, FOI_cal, aveUptake_cal, cat2uptake_cal, simParamsFolder)
%TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 'NA', 2013,0, 2.13, 0.0024, 0.0042, 2.4, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_cal235')
%%%%%%%%%%%%%%%%%%%%%%%%
% Old Function Calls:
% TBsimulation_july6(folderName, logComment,durationYrs,
%   numberPpl, plotResolution, plotUnit, recordUnit, startRecordingTime,
%   loadBurnInStr, startScenarioYr, simParamsFolder)
%
% old function call;
% function TBsimulation_july6(folderName, logComment,durationYrs, numberPpl, plotResolution, plotUnit, recordUnit, startRecordingTime, loadBurnInStr, startScenarioYr,startScenarioYr2, simParamsFolder )
% TBsimulation_july6('.', 'This is the default simulation run.',180, 480,  '-r70', 1, 60, 0, 'NA', 2013,0, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta_cal235')
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sometimes you want to run with 230 years.  In that case do not make graphs.
%
% Default settings:
% TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 1, 60, 0, 'NA', 2013, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')
%   where inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta is the folder in
%   which TBsimParamsOverwriter and TBsimParamsOverwriter_post2011 lives
%
% To write burn in:
% TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 1, 60, 0, 'writeBurnIn', 2013, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')
%
% To load burn in:
% TBsimulation_july6('.', 'This is the default simulation run.',180, 480000,  '-r70', 1, 60, 0, 'loadBurnIn', 2013, 'inf0p0023_lat2p16_aveTreat_fullCatIV_empUpta')
% If folderName is '.', the output file will be created in the same folder
% as the simParamsFolder.
%
% If plotResolution is '-r0', it does not make plots.
% plotUnit is how often to plot (1 is every month, 12 is once a year,
% etc).
%
% recordUnit is how often to record the simulation output (in months) for making
% graphs, etc.  Best to keep plotUnit = 1 if recordUnit ~=1.
%
% startRecordingTime is the time period after which to start writing records
% down.  Should normally be 0, but if want output to be Stata readable put
% it after the burn in period (i.e., 1560)
%
% loadBurnInStr is a string that will record the burn in period to a .mat if it is equal to
% 'writeBurnIn' and it will use a burn in seed (a .mat file) if equal to
% 'loadBurnIn'
%
% simParamsFolder is the string of the folder in which TBsimParams lives.
% '.' is do not change directory.
%
% startScenarioYr is between 1996 and 2045.  It is the year you want to
% turn "future projections", contained in TBsimParamsOverwriter_Post2011, scenarios on.
% startScenarioYr2 is between 1996 and 2045.  It is the year you want to
% turn "future projections", contained in TBsimParamsOverwriter_Post2019, scenarios on.
% set to 0 if there is no TBsimParamsOverwriter_Post2019 file.
%
% Notes: This is a TB simulation for the costPag effectiveness analysis of
% the Xpert rapid diagnosis system for the Indian RNTCP program.
%
% State matrix elements
% 1 | sex               | 1=male, 2=femalet
% 2 | age               | integer age
% 3 | age increment     | number of months since last bday
% 4 | smoking status    | 0=not smoke, 1=smoke
% 5 | urban/rural       | 1=urban, 2=rural
% 6 | treatment counter | number of months in treatment, if in trtmnt
% 7 | health state      | 0=healthy, 1=latentSensTB, 2=latentMDR,
%       3=activesensTB, 4=activeMDR, 5=unborn, 6=dead
% 8 | treatment         | 0=none, 1=catI, 2=catII, 3=catIII, 4=catIV
% 9 | past treat. hist. | 0=none, 1=ever and failed, 2=ever and defaulted, 3 = ever and "cured" or cured
% 10| timePerInfect     | when the individual was infected with TB
% 11| timeSmokingChange | last time individual started/stopped smoking
% 12| mdrTesting        | 0 if never tested for MDR, time Period of start MDR testing if in testing
% 13| mdrTestingStatus  | 0 = default, otherwise health status of when mdr testing started
% 14
% 15
% 16| seekTreatment     | 0 = default, otherwise 1 if willing to seek treatment in RNTCP
% 17| assignedSeekTreat | 0 = default, otherwise 1 if assigned whether willing to seek treatment in RNTCP
% 18| PtrtmtQ           | 0 = default, otherwise 1 if need to wait before being put on treatment because activated this period
% 19| PtrtmtQ_private   | 0 = default, otherwise 1 if need to wait before being put on private treatment because activated this period
% 20| PprivTrtCount     | 0 = default, otherwise number of times been in private treatment
% Supporting scripts needed:
%  TBsimParams.m        Provides data structures with parameters (TBparams.X)
%  ageBracketMaker.m    Makes age brackets for correct parameter referenced.
%  simulateEvent.m      Returns vector of logicals if event happened or not
%  smokChurnSimulateEvent simulates smoking churn
%  smokingChurnMaker    Generates the smoking churn matracies from raw data
%  smokingPrevPlotter   Plots the resultant smoking prevalence from churn
%  subsetterFunc.m      function used to find the correct population subset
%  subsetterFuncInequality.m   To find the pop subset using inequalities
%  plotByAge.m          to plot by age
%  plotByTime2.m        to make health/treatment output plots
%  plotByTime3.m        to make active TB health/treatment output plots
%  tablePrinter.m       Makes .csv outputs of graphs, in a table
%
%  recNoDecoder.m       not called,  but can be used to decode recNo
%  TBsimulationWrapper  Calls TBsimulation with specified run params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');

if isempty(folderName)
    folderName = timeStamp;
elseif strcmp(folderName, '.')
    folderName = [simParamsFolder '/' strcat(logComment,'_',timeStamp)];
else
    folderName = [folderName '/' strcat(logComment,'_',timeStamp)];
end

mkdir(folderName);
diary([folderName '/TB simulation.txt'])
diary on
startTime = cputime

%comments to text output
disp(logComment)
disp(sprintf('Folder name is', simParamsFolder))
disp(sprintf('Plot resolution is %s', plotResolution))
disp(sprintf('latToAct_cal value is %s', num2str(latToAct_cal)))
disp(sprintf('oldActSlope_cal value is %s', num2str(oldActSlope_cal)))
disp(sprintf('oldActIntercept_cal value is %s', num2str(oldActIntercept_cal)))
disp(sprintf('FOI_cal value is %s', num2str(FOI_cal)))
disp(sprintf('aveUptake_cal value is %s', num2str(aveUptake_cal)))
disp(sprintf('cat2uptake_cal value is %s', num2str(cat2uptake_cal)))

%make the burnIn seed flags
writeBurnIn = strcmp(loadBurnInStr,'writeBurnIn');
loadBurnIn = strcmp(loadBurnInStr,'loadBurnIn');
writeBurnInPostTrt = strcmp(loadBurnInStr,'writeBurnInPostTrt');
loadBurnInPostTrt = strcmp(loadBurnInStr,'loadBurnInPostTrt');
LEbuilder = strcmp(loadBurnInStr,'LEbuilder');

%%%%%initialize vars%%%%%%%%
%these used to be in the function call, but I never really change them, so
%let's make them defaults.

plotUnit = 1;
recordUnit = 60;
startRecordingTime = 0;
yrsDuration = durationYrs %should be 20
startTreatmentYr = 1996;
burnInTime = 130*12;  %number of months of burn in
burnInPostTrtTime = 149*12;  %timePeriod when 2015 starts
neverOverwrote = 1;
numPpl = numberPpl
TBparams.LEbuilderNumPpl = numberPpl;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Initialize Vars%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%load time periods
totPeriods   = 12*yrsDuration;

if loadBurnIn ~= 1  && loadBurnInPostTrt ~= 1  %not need to intialize if loading burnIn seed anyway
    %%%%%%RUN PARAM VALUES %%%%%%%%
    TBsimParams;
    totLEbuilderCohorts = size(TBparams.LEbuilderCohort,2);  %number of cohorts to do
    LEbuilderCohortDuration = TBparams.lebuilderCohortDurinYrs*12; %number of months to simulate each cohort for
    if strcmp(simParamsFolder, '.') == 0
        currentDirectory = pwd;
        cd(simParamsFolder);
        TBsimParamsOverwriter;  %want to overwrite params unique to each run with correct value
        cd(currentDirectory);
    end
    
    %%%%%PPL MATRIX, RENAMED stateMat%%%%%%%
    %columns of ppl matrix are sex, age
    Psex      = 1;  % 1 is male, 2 is female
    Page      = 2;  % age
    PageIncr  = 3;  % number of months since last bday
    Psmoking  = 4;  % 0 is not smoke, 1 is smoke
    PurbanRur = 5;  % 1 is urban, 2 is rural
    PtrtmtCounter  = 6;  % dummy for counting up months in treatment categories.
    Phealth        = 7;  % 0=healthy, 1=latentSensTB, 2=latentMDR, 3=activesensTB, 4=activeMDR, 5=unborn, 6=dead
    Ptreatment     = 8;  % current treatment is 0=none, 1=catI, 2=catII, 3=catIII, 4=catIV
    PpastTreatment = 9; % 0=none, 1=ever and failed, 2=ever and defaulted, 3 = ever and "cured" or cured
    PtimePerInfect = 10; % 0=default, otherwise month number of most recent infection.
    PtimeSmokingChange = 11; % 0 = default, otherwise month num of last smoking status change
    PmdrTesting        = 12; % 0 = default, otherwise month num of when mdr testing started
    PmdrTestingStatus  = 13; % 0 = default, otherwise health status of when mdr testing started
    PactTime           = 14; % 0 = default, otherwise time when TB activated
    PmdrEvolved        = 15; % 0 = default, otherwise 1 if mdr was evolved, not transmitted.  2 if evolved in private sector.
    PseekTrtmt         = 16; % 0 = default, otherwise 1 if willing to seek treatment in RNTCP
    PassignedSeekTrtmt = 17; % 0 = default, otherwise 1 if assigned whether willing to seek treatment in RNTCP
    PtrtmtQ            = 18; % 0 = default, otherwise 1 if need to wait before being put on treatment because activated this period
    PtrtmtQ_private    = 19; % 0 = default, otherwise 1 if need to wait before being put on private treatment because activated this period
    PprivTrtCount      = 20; % 0 = default, otherwise number of times been in private treatment
    numParams = 20; % number of parameters.  This is also hard coded in TBsimParams in the LEbuilderCohort
    % Number of different codes in Phealth
    numHealthStates = 6;
    
    %%load default population.  A matrix, one obs per person holding current state.
    numPplNeeded = floor(numPpl*((1 +TBparams.constantHisPopGrowth )^180 ) ); %floor(numPpl*((1 +TBparams.constantHisPopGrowth )^yrsDuration ) );
    fprintf('Total people by end of simulation = %d\n',numPplNeeded);
    numCurrentPpl = numPpl;
    stateMat = zeros(numPplNeeded,numParams,'int16');  %matrix of individial characteristics and health states
    stateMat(:,Phealth)   = 5; %everyone starts out as not yet born;
    stateMat(:,Psex)      = 1; %everyone starts out as male.  This only matters for making recNo, does not change analysis.
    stateMat(:,PurbanRur) = 2; %everyone starts out as rural. This only matters for making recNo, does not change analysis.
    recNo = -1 * ones(numPplNeeded,(floor((totPeriods-startRecordingTime)/recordUnit)));
    privateToRNTCP = [];  %this is initalized to nothing.
    notifications = zeros(numPplNeeded, 1);  %for TBmac notifications
    privateTrtNaive = zeros(numPplNeeded, 1);
    %%%%%%%%%
    % Initialize some counter matrices
    % Counts of who dies of what in each time period
    deathsCounterMat = zeros(totPeriods, 23);  %PARTIALLY COMMENTED OUT FOR SPEED
    numPplBorn = zeros(totPeriods,1);
    detectedPpl_private = [];
    deathsPostBurnIn = zeros(totPeriods-burnInTime+1,3);  %THROWS ZEROS FOR SPEED
    % Number of people who get transmitted to in each time period
    incidenceMat = zeros(totPeriods,3);
    % Total people and people put into treatment
    diagnosisCounterMat = zeros(totPeriods, 12);
    diagnosedPool = zeros(totPeriods, 6);
    diagnosedAges = zeros(totPeriods, 21);
    diagnosisTimeMat = zeros(totPeriods, 62);
    diagnosisTimeDist = zeros(totPeriods,20);
    mdrTestPpl = zeros(totPeriods,24);
    % for private clinic
    privateCounterMat = zeros(totPeriods, 5);
    %force of infection
    forceOfInfection = zeros(floor(totPeriods),1);  %COMMENTED OUT FOR SPEED
    % Activation
    activaCounterMat = zeros(totPeriods, 8);
    activaCountSex = zeros(totPeriods, 16); %COMMENTED OUT FOR SPEED
    latentCount = zeros(totPeriods, 16); %COMMENTED OUT FOR SPEED
    activaToday = zeros(totPeriods, 2); %COMMENTED OUT FOR SPEED
    % load ds and mdr transmission (and MDR activiation) counter matracies
    MdrMethod        = zeros(totPeriods-burnInTime+1,4);
    MdrActivation    = zeros(totPeriods-burnInTime+1,1);
    diseaseCasesMat  = zeros(totPeriods-burnInTime+1,8);
    numRefPpl  = zeros(totPeriods-burnInTime+1,3);
    numQalyPpl       = zeros(19,totPeriods-burnInTime+1);
    costsMat = zeros(totPeriods-burnInTime+1,1);
    LEnumQalyPpl = zeros(19,LEbuilderCohortDuration);  %for LEbuilder
    LEcostsMat = zeros(LEbuilderCohortDuration,1);  %for LEbuilder
    LEcostsMatBig = zeros(LEbuilderCohortDuration,totLEbuilderCohorts);   %for LEbuilder
    lastStateBig = zeros(totLEbuilderCohorts, 20);  %for LEbuilder
    LEnumQalyPpl_males = zeros(19, LEbuilderCohortDuration*(totLEbuilderCohorts/2)); %for LEbuilder
    LEnumQalyPpl_females = zeros(19, LEbuilderCohortDuration*(totLEbuilderCohorts/2)); %for LEbuilder
    trtCountMat = zeros(totPeriods,2);
    kids_TBstateCount= zeros(totPeriods,7); %COMMENTED OUT FOR SPEED
    % %load death ages matrices
    nosmokeDeathAgesMat      = zeros(totPeriods-burnInTime+1, 101);   %COMMENTED OUT FOR SPEED
    smokeDeathAgesMat      = zeros(totPeriods-burnInTime+1, 101); %COMMENTED OUT FOR SPEED
    nosmokeDSTBdeathAgesMat  = zeros(totPeriods-burnInTime+1, 101); %COMMENTED OUT FOR SPEED
    smokeDSTBdeathAgesMat  = zeros(totPeriods-burnInTime+1, 101); %COMMENTED OUT FOR SPEED
    nosmokeMDRTBdeathAgesMat = zeros(totPeriods-burnInTime+1, 101); %COMMENTED OUT FOR SPEED
    smokeMDRTBdeathAgesMat = zeros(totPeriods-burnInTime+1, 101); %COMMENTED OUT FOR SPEED
    cat1MDRdead = zeros(totPeriods,4); %COMMENTED OUT FOR SPEED
    cat2MDRdead = zeros(totPeriods,4); %COMMENTED OUT FOR SPEED
    deathAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 12);  %for TBMAC
    deathTBAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 12); %first two cols are children
    aliveAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 12);
    aliveTBAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 12);
    aliveMDRAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 4);
    aliveLatAgesMat_TBmac = zeros(totPeriods-burnInTime+1, 4);
    notification_TBmac = zeros(totPeriods,7); 
    numScreened = zeros(totPeriods,4);
    %load DOTS performance measures mat
    DOTSperformanceMeas = zeros(totPeriods-burnInTime+1, 11);
    %load past treatment numbers mat
    pastTreatment = zeros(totPeriods-burnInTime+1, 5);
    catInonSuccessTypes = zeros(totPeriods-burnInTime+1, 3);
    catIInonSuccessTypes = zeros(totPeriods-burnInTime+1, 6);
    %treatment demogs
    treatmentDemog = zeros(totPeriods-burnInTime+1, 9);
    numPrivateClinics = zeros(totPeriods, 8);
    hasTBin2003 = zeros(20,20);
    ignoredTBmat = zeros(totPeriods, 1);
    %for model reduction
    %monthlyAgeHealth = zeros(totPeriods-burnInTime+1,140);
    
    if durationYrs == 230
        %                       Time period           ages sex smoking
        aliveAgesMat = zeros(totPeriods-burnInTime+1, 101,  2,   2);
        deathAgesMat = zeros(totPeriods-burnInTime+1, 101,  2,   2);
        %                       Time period                 monsSinceActivated smoking
        aliveActTBAgesMat = zeros(totPeriods-burnInTime+1, 350, 2, 2);
        aliveActTBAgesMat_after = zeros(totPeriods-burnInTime+1, 350, 2, 2);
        ageFlagMat = zeros(size(stateMat,1),1);
        MonSinceAct = zeros(size(stateMat,1),1);
        MonSinceLat = zeros(size(stateMat,1),1);
        MonSinceActTruncations = zeros(totPeriods-burnInTime+1, 1);
        cohortFlagMat = zeros(size(stateMat,1), 1);
        %average age of smokers and nonsmokers
        aveSmokerNonAge = zeros(totPeriods-burnInTime+1, 2);
    end
    monthlyActOutcomes = zeros(totPeriods-burnInTime+1,12);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%LOAD PPL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if LEbuilder == 0
        %populate personal characteristics
        display('populating personal characterisitcs')
        % First, assign an age
        stateMat(1:numPpl,Page) = floor(interp1(TBparams.cumpercAgeBrac, TBparams.cumpercAgeBracLimits, rand(numPpl,1)));
        % Find sex based on age bracket
        stateMat(1:numPpl,Psex) = simulateEvent(1-TBparams.percMale, 1, TBparams.percMaleAgeBrac, stateMat(1:numPpl,Page), []) + 1;
        % Smoking
        stateMat(1:numPpl,Psmoking) = simulateEvent(TBparams.smokingPrev, 0, TBparams.smokingPrevAgeBrac, stateMat(1:numPpl,Page), stateMat(1:numPpl,Psex));
        % Urban/rural
        stateMat(1:numPpl,PurbanRur) = ( rand(numPpl,1) > TBparams.urbanFrac(stateMat(1:numPpl,Psex))' ) + 1;
        toc
        %%%%%%%%PLOT INITIAL PEOPLE CHARS%%%%%%%%%%
        if (strcmp(plotResolution,'-r0') == 0 )
            display('plotting initial cohort')
            plotByAge(stateMat, Page, [Psex Psmoking PurbanRur]);
            print('-dpng',plotResolution,[folderName '/startingDemogsPlot.png']);
            close all
        end
        
        %%%%%%%%%%%%%%%%%%%%%Initialize disease prevalence%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Initializing TB prevalence');
        % Start with everyone healthy
        stateMat(1:numPpl,Phealth) = 0;
        % latent sensTB population
        latentTB = simulateEvent(TBparams.latentTBprev,0,TBparams.latentTBageBrac,stateMat(1:numPpl,Page),[]);
        stateMat(latentTB,Phealth) = 1;
        stateMat(latentTB,PtimePerInfect) = 1; %set infection time as 1, so not 0 (which is default)
        %%%%%% active TB population
        % Nonsmoking
        nonSmoking = find(stateMat(1:numPpl,Psmoking) == 0); % Indices
        activeTB = simulateEvent(TBparams.ruralTBprev, 0, TBparams.ruralTBprevAgeBrac, stateMat(nonSmoking,Page), stateMat(nonSmoking,Psex)); % Logical
        stateMat(nonSmoking(activeTB),Phealth) = 3;
        stateMat(nonSmoking(activeTB),PtimePerInfect) = 1; %set infection time as 1, so not 0 (which is default)
        
        % smoking
        smoking = find(stateMat(1:numPpl,Psmoking) == 1);
        activeTB = simulateEvent(TBparams.ruralTBprev * TBparams.TBprevSmokingMultiplier, 0, TBparams.ruralTBprevAgeBrac, stateMat(smoking,Page), stateMat(smoking,Psex));
        stateMat(smoking(activeTB),Phealth) = 3;
        stateMat(smoking(activeTB),PtimePerInfect) = 1; %set infection time as 1, so not 0 (which is default)
        disp('Finished initializing TB prevalence')
        toc
        %%%%%%%%PLOT INITIAL HEALTH CHARS%%%%%%%%%%
        if (strcmp(plotResolution,'-r0') == 0 )
            inithealthMat(1:numPpl,numHealthStates+2) = stateMat(1:numPpl,Page);
            for i = 0:numHealthStates
                inithealthMat(1:numPpl,i+1) = (stateMat(1:numPpl,Phealth) == i);
            end
            plotByAge(inithealthMat, numHealthStates+2, [1 2 3 4 5 7]);
            title('Initial Health Distribution');
            print('-dpng',plotResolution,[folderName '/startingPrevalencePlot.png']);
            close all
            
            inithealthMat(1:numPpl,2) = inithealthMat(1:numPpl,2) + inithealthMat(1:numPpl,3);
            inithealthMat(1:numPpl,4) = inithealthMat(1:numPpl,4) + inithealthMat(1:numPpl,5);
            inithealthMat(1:numPpl,6) = inithealthMat(1:numPpl,7);
            plotByAge(inithealthMat, numHealthStates+2, [1 2 4 6]);
            title('Initial Health Distribution (small)');
            print('-dpng',plotResolution,[folderName '/startingPrevalencePlotSmall.png']);
            close all
        end
    end %end LEbuilderif
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%TIME LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeLoopTimer = zeros(14,1);
timeLoopTimerName = cell(14,1);
post2011RunYet = 0;
post2019RunYet = 0;
TBparams.numBirths = dlmread('numPplBorn.csv', ',' ,1,0);  %for hard coding the number of births

%time loop.  Loop for totPeriods in increments of 1 month
timePeriod = 1;
if LEbuilder ==1
    cohortNum = 0;
    origTotPeriods = totPeriods;
    timePeriod = origTotPeriods+1;
    totPeriods = origTotPeriods + LEbuilderCohortDuration;  %timePeriod + 1200 + 2;
end

while timePeriod <= totPeriods
    timeLoopIndex = 1;
    if writeBurnIn == 1
        %save all the burn in stuff to make a burnIn seed
        if timePeriod == burnInTime - 1
            oldPeriods = totPeriods;
            clear logComment folderName totPeriods durationYrs simParamsFolder plotResolution plotUnit startRecordingTime startScenarioYr startScenarioYr2 writeBurnIn loadBurnIn aveUptake_cal cat2uptake_cal
            save('burnInTester', '-v7.3');
            plotResolution = '-r0';
            plotUnit = 1;
            startRecordingTime = 0;
            writeBurnIn = 1;
            loadBurnIn = 0;
            totPeriods = oldPeriods;
            break
        end
    end
    
    if writeBurnInPostTrt == 1  %this is writing a mat file in 2015
        %save all the stuff till 2014 to make a burnIn seed
        if timePeriod == burnInPostTrtTime - 1
            oldPeriods = totPeriods;
            clear logComment folderName totPeriods durationYrs simParamsFolder plotResolution plotUnit startRecordingTime startScenarioYr startScenarioYr2 writeBurnIn loadBurnIn aveUptake_cal cat2uptake_cal
            save('burnInTesterPostTrt', '-v7.3');
            plotResolution = '-r0';
            plotUnit = 1;
            startRecordingTime = 0;
            writeBurnIn = 1;
            loadBurnIn = 0;
            totPeriods = oldPeriods;
            break
        end
    end
    
    %load the burn in stuff if it's time to do that and user specified it
    if (timePeriod == 1 && loadBurnIn == 1)
        fprintf('loading burnIn seed\n');
        clearvars -except logComment folderName durationYrs totPeriods simParamsFolder plotResolution plotUnit startRecordingTime startScenarioYr startScenarioYr2 writeBurnIn loadBurnIn aveUptake_cal cat2uptake_cal
        load burnInTester;
        TBparams.numBirths = dlmread('numPplBorn.csv', ',' ,1,0);  %for hard coding the number of births
        privateToRNTCP = [];  %this is initalized to nothing.
        timePeriod = burnInTime - 1;
    end
    
    if (timePeriod == 1 && loadBurnInPostTrt  == 1)
        fprintf('loading post-treatment burnIn seed\n');
        clearvars -except logComment folderName durationYrs totPeriods simParamsFolder plotResolution plotUnit startRecordingTime startScenarioYr startScenarioYr2 writeBurnIn loadBurnIn aveUptake_cal cat2uptake_cal
        load burnInTesterPostTrt;
        TBparams.numBirths = dlmread('numPplBorn.csv', ',' ,1,0);  %for hard coding the number of births
        timePeriod = burnInPostTrtTime - 1;
    end
    
    if (neverOverwrote == 1 && (loadBurnIn == 1 || loadBurnInPostTrt  == 1 || LEbuilder == 1))
        %comments to text output
        disp('loading TBsimParamsOverwriter')
        disp(sprintf('latToAct_cal value is %s', num2str(latToAct_cal)))
        disp(sprintf('FOI_cal value is %s', num2str(FOI_cal)))
        disp(sprintf('aveUptake_cal value is %s', num2str(aveUptake_cal)))
        disp(sprintf('cat2uptake_cal value is %s', num2str(cat2uptake_cal)))
        
        %reload the overwriter
        if strcmp(simParamsFolder, '.') == 0
            currentDirectory = pwd;
            cd(simParamsFolder);
            TBsimParamsOverwriter;  %want to overwrite params unique to each run with correct value
            cd(currentDirectory);
            
            if strncmp(logComment,'r',1)  %run sensitvity scenario
                disp('Using a sensitivity analysis scenario:')
                sensitivityScenarioNum = str2num(logComment(2:3))
                sensitivityScenario
            end
        end
        neverOverwrote = 0;
    end
    
    if LEbuilder == 0
        if (mod(timePeriod,10)==0)
            fprintf('Time period %i\n',timePeriod);
        else
            fprintf('.');
        end
    end
    %if post 2011 we want to turn particular scenarios on.  changed to 2013.
    if startScenarioYr ~= 0
        if (timePeriod >= TBparams.yearToMon(find(TBparams.yearToMon(:,1)==startScenarioYr),2) && post2011RunYet ==0)
            fprintf('\n Changing Overwriter to Post one.  Time Period: %i (calendar year) %i \n', timePeriod, TBparams.yearToMon(find(TBparams.yearToMon(:,2)==timePeriod),1));
            if strcmp(simParamsFolder, '.') == 0                
                currentDirectory = pwd;
                cd(simParamsFolder);
                TBsimParamsOverwriter_post2011;  %want to overwrite params unique to each run with correct value
                cd(currentDirectory);
                
                if strncmp(logComment,'r',1)  %reload the sensitivity scenario
                    disp('Using a sensitivity analysis scenario:')
                    %run sensitvity scenario
                    sensitivityScenarioNum = str2num(logComment(2:3))
                    sensitivityScenario
                end
            end
            post2011RunYet = 1;
        end
    end
    
    if startScenarioYr2 ~= 0
        if (timePeriod >= TBparams.yearToMon(find(TBparams.yearToMon(:,1)==startScenarioYr2),2) && post2019RunYet ==0)
            fprintf('\n Changing Overwriter (2nd time).  Time Period: %i (calendar year) %i \n', timePeriod, TBparams.yearToMon(find(TBparams.yearToMon(:,2)==timePeriod),1));
            if strcmp(simParamsFolder, '.') == 0
                currentDirectory = pwd;
                cd(simParamsFolder);
                TBsimParamsOverwriter_post2019;  %want to overwrite params unique to each run with correct value
                cd(currentDirectory);
                %DEBUGGING
                TBparams.DotsPlusRampUpPeriod
                TBparams.DotsPlusFullCoverage
                TBparams.DotsPlusRampUpSequence
            end
            post2019RunYet = 1;
        end
    end
    %in 1996 we want to seed some of the population with MDR TB
    if timePeriod == burnInTime && LEbuilder == 0
        DSTB_forMDR_lat = subsetterFunc(stateMat,[-1],Phealth,[1]); % convert some of the lat DS TB to lat MDR TB
        hasMDR_log_lat = rand(length(DSTB_forMDR_lat),1) < TBparams.MDRseedIn1996_lat ;
        stateMat(DSTB_forMDR_lat(hasMDR_log_lat), Phealth) = stateMat(DSTB_forMDR_lat(hasMDR_log_lat), Phealth) + 1;  %change them to MDR instead of DS
        DSTB_forMDR_act = subsetterFunc(stateMat,[-1],Phealth,[3]); % convert some of the act DS TB to act MDR TB
        hasMDR_log_act = rand(length(DSTB_forMDR_act),1) < TBparams.MDRseedIn1996_act ;
        stateMat(DSTB_forMDR_act(hasMDR_log_act), Phealth) = stateMat(DSTB_forMDR_act(hasMDR_log_act), Phealth) + 1;  %change them to MDR instead of DS
    end
    
    timerStart = cputime;
    
    
    %DEBUGGING
    % sensSpecDebugMat = [
    %     TBparams.SSsensit,TBparams.SSspec,TBparams.DSTsensit, TBparams.DSTspec;
    %     TBparams.geneXsensDS,TBparams.geneXspecDS,TBparams.geneXsensMDR,TBparams.geneXspecMDR
    % ]
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%BIRTHS%%%%%%%%%%%%%%%%%%%%%%%%
    
    % have births = number of deaths + (  TBparams.popGrowthPerc/(3650*12) * lastMonPop  )
    births = find(stateMat(:,Phealth) == 6);
    if LEbuilder == 0
        %get newly born people to make population grow
        if (timePeriod > burnInTime+1 && numberPpl == 200000)
            lastNewlyBornPerson =(TBparams.numBirths(timePeriod,1) - size(births,1))+numCurrentPpl;  %this is same number of births for all interventions
            if (lastNewlyBornPerson < numCurrentPpl)  %error if neg births
                disp('error: negative pop growth.  Number of births supposed to have, num dead ppl to make up for, num current people:')
                TBparams.numBirths(timePeriod,1)
                size(births,1)
                numCurrentPpl
                error('Error: negative pop growth')
            end            
        else
            %first find the appropriate pop growth for this period
            annualPopGrowthRate = interp1( TBparams.hisPopGrowthYrs, TBparams.hisPopGrowthPerc, startTreatmentYr + subplus((timePeriod - burnInTime)/12) );
            lastNewlyBornPerson = floor(   numCurrentPpl * ((1 + annualPopGrowthRate)^(1/12) ) ); %this is the normal birth process
        end
        
        births = [births; (numCurrentPpl+1:lastNewlyBornPerson)']; % Indices
        numPplBorn(timePeriod,1) = size(births,1);  %number of births
        %change numCurrentPpl to the right number
        numCurrentPpl = lastNewlyBornPerson;
        
        % Reset their lives (except sex and urban/rural)
        stateMat(births,1:numParams) = 0;  %shorthand, makes 0 for all params for all born (dimensions not look right but okay for matlab)
        %sex and urban/rural need to be reassigned randomly
        stateMat(births,Psex)      = simulateEvent(1-TBparams.percMale, 1, TBparams.percMaleAgeBrac, stateMat(births,Page), []) + 1;
        stateMat(births,PurbanRur) = ( rand(length(births),1) > TBparams.urbanFrac(stateMat(births,Psex))' ) + 1;
        notifications(births, 1) = 0;  %reset notifications for TBmac
        privateTrtNaive(births, 1) = 0;  %reset private trt naiveness for TBmac
 
        if mod(timePeriod,12)==0
            notifications(:,1) = 0;  %reset all notification counters for tbmac
        end
        
        %debug and calibration: make a cohort made up of the people born in three years period after burn in
        if yrsDuration ==230
            if (timePeriod >= burnInTime+1  && timePeriod <= burnInTime +36 )
                cohortFlagMat(births,1) = 1;
            else
                cohortFlagMat(births,1) = 0;
            end
        end
        % end debug and calibration
    else  %LEbuilder is 1
        stateMat(births,1:numParams) = 0;  %if LEbuilder is 1, reset everything but keep them dead
        stateMat(births,1:Phealth) = 6;
        if timePeriod == origTotPeriods+1 %LEbuilder is 1 and need to put in a new cohort (put in a new one every 100 yrs)
            cohortNum = cohortNum+1;
            
            stateMat = repmat(TBparams.LEbuilderCohort{cohortNum},TBparams.LEbuilderNumPpl,1);
            %stateMat = TBparams.LEbuilderCohort{cohortNum};
            
            %reset costs and qalys
            costsMat = zeros(totPeriods-burnInTime+1,1);
            LEnumQalyPpl = zeros(19,LEbuilderCohortDuration); %zeros(19,1200);
            LEcostsMat = zeros(LEbuilderCohortDuration,1);  %zeros(1200,1);
            numCurrentPpl = TBparams.LEbuilderNumPpl;  %instead of 10000
        end
    end
    
    if timePeriod >= burnInTime
        aliveAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(stateMat(1:numCurrentPpl, Page),[0;10;15;20;30;40;50;60;70;80;90;110])';
        hasTB = subsetterFunc(stateMat,[-1],Phealth,[3,4]);
        aliveTBAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(stateMat(hasTB, Page),[0;10;15;20;30;40;50;60;70;80;90;110])';
        hasMDR = subsetterFunc(stateMat,[-1],Phealth,[4]);
        aliveMDRAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(stateMat(hasMDR, Page),[0;10;15;110])';
        haslatent = subsetterFunc(stateMat,[-1],Phealth,[1,2]);
        aliveLatAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(stateMat(haslatent, Page),[0;10;15;110])';
    end
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'births';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%%%%INCREMENT AGE %%%%%%%%%%%%%%%%%%
    %age = age + 1 month
    stateMat(1:numCurrentPpl,PageIncr) = stateMat(1:numCurrentPpl,PageIncr) + 1;
    agedAyear = stateMat(1:numCurrentPpl,PageIncr) == 12; % logical
    stateMat(agedAyear,Page) = stateMat(agedAyear,Page) + 1;
    stateMat(agedAyear,PageIncr) = 0;
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'aging';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%%%%%%%%%%%%%SMOKERS%%%%%%%%%%%%%%%%%%%%%%%%
    %smoking status can change yearly (not monthly)
    %     if mod(timePeriod,12)== 0
    %         %make groups
    %         smokers    = subsetterFunc(stateMat,[-1],Psmoking,[1]); % indices
    %         nonsmokers = subsetterFunc(stateMat,[-1],Psmoking,[0]); % indices
    %         urbanSmokers    = subsetterFunc(stateMat,   smokers,PurbanRur,[1]); % these are all indices, too
    %         urbanNonSmokers = subsetterFunc(stateMat,nonsmokers,PurbanRur,[1]);
    %         ruralSmokers    = subsetterFunc(stateMat,   smokers,PurbanRur,[2]);
    %         ruralNonSmokers = subsetterFunc(stateMat,nonsmokers,PurbanRur,[2]);
    %         % get people to start or stop smoking, depending on their groups
    %         churnSubset = {urbanSmokers, urbanNonSmokers, ruralSmokers, ruralNonSmokers};
    %         churnProbs  = {TBparams.oneMonQuittingUrbanProb,  TBparams.oneMonStartingUrbanProb, TBparams.oneMonQuittingRuralProb, TBparams.oneMonStartingRuralProb};
    %         churnFactor = {TBparams.smoker_quit_ratio , TBparams.nonsmoker_start_ratio,TBparams.smoker_quit_ratio, TBparams.nonsmoker_start_ratio};
    %         smokingChangeResult = {0, 1, 0, 1};
    %         for i = 1:4
    %             %four categories of smoking/nonsmoking: smoking 0 yrs, 1 yrs, 2 yrs, 3yrs, nonsmoking etc.
    %             yearsSinceSmoked = floor(double(timePeriod - stateMat(churnSubset{i},PtimeSmokingChange))/12);
    %             yearsSinceSmoked(yearsSinceSmoked > 4) = 4;
    %             changeStat = smokChurnSimulateEvent(churnProbs{i}, churnFactor{i}, yearsSinceSmoked, 0, TBparams.smokingChurnAgeBrac, stateMat(churnSubset{i},Page), stateMat(churnSubset{i},Psex));
    %             stateMat(churnSubset{i}(changeStat),Psmoking) = smokingChangeResult{i}; %quit smoking
    %             stateMat(churnSubset{i}(changeStat),PtimeSmokingChange) = -2; % flag that their timesSmokingChange should be changed at the end
    %         end
    %         % Update stateMat(:,PtimeSmokingChange)
    %         % People who have a new status need to have their timeSmokingChange
    %         newStatus       =  stateMat(:,PtimeSmokingChange) == -2; % logical
    %         % Also some nonsmokers under 18 have memory problems apparently
    %         % (they are called retained)
    %         under18notSmoke = (stateMat(:,Page) < 18 & stateMat(:,Psmoking) == 0);  % logical
    %         probRetained    = double((18.0 - stateMat(under18notSmoke,Page)))/18.0; % prob of not incremented according to this func
    %         retainedUnder18notSmoke = rand(sum(under18notSmoke),1) < probRetained ; % simulate whether retained
    %         retained = under18notSmoke;  %logical, length of stateMat
    %         retained(under18notSmoke) = retainedUnder18notSmoke; %length of statemat and correct people
    %         % Update the timeSmokingChange
    %         stateMat(retained, PtimeSmokingChange) = stateMat(retained,PtimeSmokingChange) + 12; % For retained people, increment by a year
    %         stateMat(newStatus,PtimeSmokingChange) = timePeriod; % For the new status people, set smoking change time to current time period
    %     end
    %
    %     %FOR DEBUGGING ONLY TO LOOK AT PREVALENCE
    %     if (mod(timePeriod,240)==0 && timePeriod >= burnInTime)
    %         if (strcmp(plotResolution,'-r0') == 0 )
    %             smokingPrevPlotter(stateMat, numCurrentPpl, timePeriod, Page, Psex, Psmoking, plotResolution, folderName, 0);
    %             smokingPrevPlotter(stateMat, numCurrentPpl, timePeriod, Page, Psex, Psmoking, plotResolution, folderName, 1);
    %         end
    %     end
    %
    %    if durationYrs == 230
    %       if timePeriod >= burnInTime
    %          aveSmokerNonAge(timePeriod-burnInTime+1, 1) = mean(stateMat(find(stateMat(:, Psmoking) == 1), Page));
    %          aveSmokerNonAge(timePeriod-burnInTime+1, 2) = mean(stateMat(find(stateMat(:, Psmoking) == 0), Page));
    %          aveSmokerNonAge(timePeriod-burnInTime+1, 3) = mean(stateMat(find(stateMat(:, Psmoking) == 1 & stateMat(:,Page) >= 30 ), Page) );
    %          aveSmokerNonAge(timePeriod-burnInTime+1, 4) = mean(stateMat(find(stateMat(:, Psmoking) == 0 & stateMat(:,Page) >= 30 ), Page));
    %       end
    %    end
    %
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'smoking';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    
    %%%%%%%%%%%%%%%%DEATH%%%%%%%%%%%%%%%%%%%%%%
    %find the appropriate mortality numbers to use;
    deathYr = floor((timePeriod-burnInTime)/12) + 7;
    if deathYr < 7,  deathYr = 7;  end % Use the 1996 numbers for everything pre-1996 (deathYr = 1 means 1990 mortality)
    if deathYr > 20, deathYr = 20; end % Use the 2009 numbers for everything post-2009
    % natural death rates requires change to 1990 or 2000 death year (deathYr = 1 or 11)
    if TBparams.fixedDeathYr ~= 0
        deathYr = TBparams.fixedDeathYr;  %this is for debug or calibration only
    end
    
    deathAges  = []; deathAgesTB = []; %initialize for TB MAC
    
    %for death rates debug
    % Count the number of people alive
    if durationYrs == 230
        if timePeriod >= burnInTime
            for fish = 1:numCurrentPpl
                if cohortFlagMat(fish,1) == 1
                    aliveAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) ...
                        = aliveAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) + 1;
                    if stateMat(fish,Phealth) == 1 || stateMat(fish, Phealth) == 2  %count latent people in the cohort
                        MonSinceLat(fish) = MonSinceLat(fish) + 1;
                    end
                    if stateMat(fish, Phealth) == 3 || stateMat(fish, Phealth) == 4  %count them in the actTBmatrix
                        MonSinceAct(fish) = MonSinceAct(fish) + 1;
                        if MonSinceAct(fish) > size(aliveActTBAgesMat,2)
                            MonSinceActTruncations(timePeriod-burnInTime+1, 1) = MonSinceActTruncations(timePeriod-burnInTime+1, 1)+1;
                            MonSinceAct(fish) = size(aliveActTBAgesMat,2);
                        end
                        if stateMat(fish,Page) >= 30
                            ageFlagMat(fish) = 1;
                        end
                        aliveActTBAgesMat(timePeriod-burnInTime+1, MonSinceAct(fish), stateMat(fish,Psmoking)+1, ageFlagMat(fish)+1) ...
                            = aliveActTBAgesMat(timePeriod-burnInTime+1, MonSinceAct(fish), stateMat(fish,Psmoking)+1, ageFlagMat(fish)+1) + 1;
                    end
                end
            end
        end
    end
    % end for natural death rates debug
    
    % ---------- NATURAL CAUSES --------------------------------------------------
    %dying from natural causes (age and gender dependent)
    natDeathnonActTB = subsetterFunc(stateMat,[-1],Phealth,[0,1,2]);
    %if on perfect treatment need to suffer from natural causes of death
    natDeathCatI = []; natDeathCatII=[];
    if any(TBparams.CatIdeath == 0 )
        natDeathCatI = subsetterFunc(stateMat,[-1],Ptreatment,[1]);
    end
    if any(TBparams.CatIIdeath == 0 )
        natDeathCatII = subsetterFunc(stateMat,[-1],Ptreatment,[2]);
    end
    if ( isempty(natDeathCatI) && isempty(natDeathCatII) )
        nonActTB = natDeathnonActTB;
    else
        nonActTB = sort(unique([natDeathnonActTB;natDeathCatI; natDeathCatII]));
    end
    
    %nonSmoking
    %    nonSmoking = subsetterFunc(stateMat, nonActTB, Psmoking,[0]); % Nonsmoking
    nonSmoking =  nonActTB;   %not doing smoking anymore, just do everyone
    dead = simulateEvent(TBparams.mortNonsmoking{deathYr}, 0, TBparams.baseMortAgeBrac, stateMat(nonSmoking,Page), stateMat(nonSmoking,Psex));
    deathsCounterMat(timePeriod, 1) = sum(dead);  %count their death
    stateMat(nonSmoking(dead),Phealth) = 6; % dies
    deathAges = [deathAges ;stateMat(nonSmoking(dead),Page)];
    if timePeriod >= burnInTime
        %initialize runningDead vector
        runningDead    = zeros(size(stateMat,1),1);
        runningDSDead  = zeros(size(stateMat,1),1);
        runningMDRDead = zeros(size(stateMat,1),1);
        runningDead(nonSmoking(dead)) = 1;
    end
    % BEGIN for natural death rates debug
    % Count the number of people died (these will all be nonsmokers)
    if durationYrs == 230
        if timePeriod >= burnInTime
            for fish = nonSmoking(dead)'
                if cohortFlagMat(fish,1) == 1
                    deathAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) ...
                        = deathAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) + 1;
                end
            end
        end
    end
    % END for natural death rates debug
    
    %smoking  NOT NEED SINCE NOT DOING SMOKING ANYMORE
    %     smoking = subsetterFunc(stateMat, nonActTB, Psmoking,[1]); % Smoking
    %     dead = simulateEvent(TBparams.mortSmoking{deathYr}, 0, TBparams.baseMortAgeBrac, stateMat(smoking,Page), stateMat(smoking,Psex));
    % %     deathsCounterMat(timePeriod, 2) = sum(dead); %count their death
    %     stateMat(smoking(dead),Phealth) = 6; % dies
    %     if timePeriod >= burnInTime
    %         runningDead(smoking(dead)) = 1;
    %     end
    %     % BEGIN for natural death rates debug
    %     % Count the number of people died (these will all be smokers)
    %     if durationYrs == 230
    %         if timePeriod >= burnInTime
    %             for fish = smoking(dead)'
    %                 if cohortFlagMat(fish,1) == 1
    %                     deathAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) ...
    %                         = deathAgesMat(timePeriod-burnInTime+1, stateMat(fish,Page)+1, stateMat(fish,Psex), stateMat(fish,Psmoking)+1) + 1;
    %                 end
    %             end
    %         end
    %     end  END NOT DOING SMOKING ANYMORE
    % END for natural death rates debug
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'natural deaths';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    % ---------- TB-RELATED --------------------------------------------------
    %Now with TB
    TBpeople =  find( stateMat(1:numCurrentPpl,Phealth) == 3 | stateMat(1:numCurrentPpl,Phealth) == 4 ) ;  %active sens/MDR TB
    MDRpeople =  find( stateMat(1:numCurrentPpl,Phealth) == 4 ) ;  %active MDR TB
    DSpeople = find( stateMat(1:numCurrentPpl,Phealth) == 3 ) ;  %active DS TB

    %for people who have untreated TB
    TBnoTrtmt = subsetterFunc(stateMat,TBpeople,Ptreatment,0);    %active sens/MDR TB people with no treatment
    MdrCatI = subsetterFunc(stateMat,MDRpeople,Ptreatment,1);    %active MDR TB people with CatI treatment
    MdrCatII = subsetterFunc(stateMat,MDRpeople,Ptreatment,2);    %active MDR TB people with CatII treatment
    TBnoTrtmt = sort([TBnoTrtmt;MdrCatI;MdrCatII]);
    %...who do not smoke
    %    nonsmokTBnoTrtmt = subsetterFunc(stateMat,TBnoTrtmt,Psmoking,[0]);
    nonsmokTBnoTrtmt = TBnoTrtmt;  % NOT DOING SMOKING ANYMORE
    dead = simulateEvent(TBparams.nonDOTsTBmortNonSmok{deathYr}, 0, TBparams.baseMortAgeBrac, stateMat(nonsmokTBnoTrtmt,Page), stateMat(nonsmokTBnoTrtmt,Psex));
    
    %%TB MAC ONLY: DEATH ADJUSTOR.  Half of these slated for death do not die, actually self cure.
    actuallyNotdead = rand(length(nonsmokTBnoTrtmt),1) < 0.43;
    selfCures = (dead == 1 & actuallyNotdead == 1);
    dead = (dead == 1 & actuallyNotdead == 0);
    stateMat(nonsmokTBnoTrtmt(selfCures),Phealth) = 0;  %self cures!
    timeSinceHadTB = timePeriod - stateMat(nonsmokTBnoTrtmt(selfCures), PactTime);  %count these people
    ignoredTBtime = [timeSinceHadTB];  %write down duration of ignored TB
    %end the TB mac death adjustor  
    
    deathsCounterMat(timePeriod,3) = sum(dead); %count their death
    timeSinceHadTB = timePeriod - stateMat(nonsmokTBnoTrtmt(dead), PactTime);  %count these people
    ignoredTBtime = [ignoredTBtime;timeSinceHadTB];  %write down duration of ignored TB
    deathsCounterMat(timePeriod,13) = size( subsetterFunc(stateMat,nonsmokTBnoTrtmt(dead),Phealth,[4]) ,1); %count their death
    deathsCounterMat(timePeriod, 19) = size(intersect(nonsmokTBnoTrtmt(dead),detectedPpl_private),1);%count their death
    deathsCounterMat(timePeriod,20) = size( subsetterFunc(stateMat,nonsmokTBnoTrtmt(dead),PpastTreatment,[1,2]) ,1); %count their death
    deathsCounterMat(timePeriod,21) = size( subsetterFunc(stateMat,nonsmokTBnoTrtmt(dead),PpastTreatment,[1,2,3]) ,1); %count their death
    deathsCounterMat(timePeriod,22) = size( subsetterFuncInequality(stateMat,nonsmokTBnoTrtmt(dead),PprivTrtCount, 1, 20) ,1); %count their death
    
%     if (TBMAC1_increaseTrt == 1 ||  TBMAC2_initDefReduction == 1||  TBMAC3_XpertReplacesSmear == 1 || TBMAC4_activeCaseFinding==1 || TBMAC5_preventativeTherapyForLat ==1 )
%         %hypothetical deaths from non-TB causes only
%         hypothetical_dead = simulateEvent(TBparams.mortNonsmoking{deathYr}, 0, TBparams.baseMortAgeBrac, stateMat(TBpeople,Page), stateMat(TBpeople,Psex));
%         deathsCounterMat(timePeriod, 23) = sum(hypothetical_dead);
%     end
%     
    stateMat(nonsmokTBnoTrtmt(dead),Phealth) = 6;  %dies
    deathAgesTB  = [deathAgesTB ;stateMat(nonsmokTBnoTrtmt(dead),Page)];
    %    trtCountMat(timePeriod, 1) = sum(dead);  %record the number of active TB, untreated deaths
    if timePeriod >= burnInTime
        runningDead(nonsmokTBnoTrtmt(dead)) = 1;
        runningDSDead(nonsmokTBnoTrtmt(dead)) = 1;
    end
    
    %     cat1MDRdead(timePeriod,1) = sum(stateMat(nonsmokTBnoTrtmt(dead),Ptreatment) == 1);
    %     cat2MDRdead(timePeriod,1) = sum(stateMat(nonsmokTBnoTrtmt(dead),Ptreatment) == 2);
    %...who smoke  NOT DOING SMOKING ANYMORE
    %     smokTBnoTrtmt = subsetterFunc(stateMat,TBnoTrtmt,Psmoking,[1]);
    %     dead = simulateEvent(TBparams.nonDOTsTBmortSmoking{deathYr}, 0, TBparams.baseMortAgeBrac, stateMat(smokTBnoTrtmt,Page), stateMat(smokTBnoTrtmt,Psex));
    % %     deathsCounterMat(timePeriod,14) = size( subsetterFunc(stateMat,smokTBnoTrtmt(dead),Phealth,[4]) ,1); %count their death
    % %     deathsCounterMat(timePeriod,4) = sum(dead);
    %     stateMat(smokTBnoTrtmt(dead),Phealth) = 6;
    %     if timePeriod >= burnInTime
    %         runningDead(smokTBnoTrtmt(dead)) = 1;
    %         runningDSDead(smokTBnoTrtmt(dead)) = 1;
    %     end  END NOT DOING SMOKING ANYMORE
    %     cat1MDRdead(timePeriod,1) = cat1MDRdead(timePeriod,1) + sum(stateMat(smokTBnoTrtmt(dead),Ptreatment) == 1);
    %     cat1MDRdead(timePeriod,4) = size(MdrCatI,1);
    %     cat2MDRdead(timePeriod,1) = cat2MDRdead(timePeriod,1) + sum(stateMat(smokTBnoTrtmt(dead),Ptreatment) == 2);
    %     cat2MDRdead(timePeriod,4) = size(MdrCatII,1);
    %death prob if in catI treatment
    catITBppl = subsetterFunc(stateMat,DSpeople,Ptreatment,[1]);  %active sens/MDR TB people with CatI treatment
    
    if ( TBMAC2b_trtSuccess == 1 && timePeriod > 1800)  %TBMAC2b
        dead = simulateEvent(TBparams.CatIdeath  *cat12DeathInitScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIageBrac, stateMat(catITBppl,Page), stateMat(catITBppl,Psex));
    else
        dead = simulateEvent(TBparams.CatIdeath, 0, TBparams.TBCatIageBrac, stateMat(catITBppl,Page), stateMat(catITBppl,Psex));
    end
    
    %%TB MAC ONLY: DEATH ADJUSTOR.  Half of these slated for death do not die, actually self cure.
    actuallyNotdead = rand(length(catITBppl),1) < 0.43;
    selfCures = (dead == 1 & actuallyNotdead == 1);
    dead = (dead == 1 & actuallyNotdead == 0);
    stateMat(catITBppl(selfCures),Phealth) = 0;  %self cures!
    stateMat(catITBppl(selfCures),Ptreatment) = 0;  %self cures and leaves trt!
    %end the TB mac death adjustor  

    %    nonSmokingDead = (stateMat(catITBppl(dead),Psmoking) == 0);  %NOT DOING SMOKING ANYMORE
    %    deathsCounterMat(timePeriod,5) = sum(nonSmokingDead);
    deathsCounterMat(timePeriod,6) = sum(dead) - deathsCounterMat(timePeriod,5);  %smoking dead
    deathsCounterMat(timePeriod,15) = size( subsetterFunc(stateMat,catITBppl(dead),Phealth,[4]) ,1); %count their death
    stateMat(catITBppl(dead),Phealth) = 6;
    deathAgesTB = [deathAgesTB ;stateMat(catITBppl(dead),Page)];
    if timePeriod >= burnInTime
        runningDead(catITBppl(dead)) = 1;
        runningDSDead(catITBppl(dead)) = 1;
    end
    %     cat1MDRdead(timePeriod,2) = sum(stateMat(catITBppl(dead),Ptreatment) == 1);
    %     cat1MDRdead(timePeriod,3 ) = size(find( stateMat(1:numCurrentPpl,Ptreatment) == 1 ),1) ;
    %death prob if in catII treatment
    catIITBppl = subsetterFunc(stateMat,DSpeople,Ptreatment,[2]);  %active sens/MDR TB people with CatII treatment
    
    if (TBMAC2b_trtSuccess == 1 && timePeriod > 1800 ) %TBMAC2b
        dead = simulateEvent(TBparams.CatIIdeath*cat12DeathInitScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIIageBrac, stateMat(catIITBppl,Page), stateMat(catIITBppl,Psex));
    else
        dead = simulateEvent(TBparams.CatIIdeath, 0, TBparams.TBCatIIageBrac, stateMat(catIITBppl,Page), stateMat(catIITBppl,Psex));
    end
    
    %%TB MAC ONLY: DEATH ADJUSTOR.  Half of these slated for death do not die, actually self cure.
    actuallyNotdead = rand(length(catIITBppl),1) < 0.43;
    selfCures = (dead == 1 & actuallyNotdead == 1);
    dead = (dead == 1 & actuallyNotdead == 0);
    stateMat(catIITBppl(selfCures),Phealth) = 0;  %self cures!
    stateMat(catIITBppl(selfCures),Ptreatment) = 0;  %self cures and leaves trt!
    %end the TB mac death adjustor  

    %    nonSmokingDead = (stateMat(catIITBppl(dead),Psmoking) == 0); %NOT DOING SMOKING ANYMORE
    %     deathsCounterMat(timePeriod,7) = sum(nonSmokingDead);
    deathsCounterMat(timePeriod,8) = sum(dead) - deathsCounterMat(timePeriod,7);  %smoking dead
    deathsCounterMat(timePeriod,16) = size( subsetterFunc(stateMat,catIITBppl(dead),Phealth,[4]) ,1); %count their death
    stateMat(catIITBppl(dead),Phealth) = 6;
    deathAgesTB = [deathAgesTB ;stateMat(catIITBppl(dead),Page)];
    if timePeriod >= burnInTime
        runningDead(catIITBppl(dead)) = 1;
        runningDSDead(catIITBppl(dead)) = 1;
    end
    %     cat2MDRdead(timePeriod,2) = sum(stateMat(catIITBppl(dead),Ptreatment) == 2);
    %     cat2MDRdead(timePeriod,3 ) = size(find( stateMat(1:numCurrentPpl,Ptreatment) == 2 ),1) ;
    
    %death prob if in catIV treatment
    catIVTBppl = subsetterFunc(stateMat,TBpeople,Ptreatment,[4]); %active TB people with CatIV treatment
    
    if (TBMAC2c_MDRtrtSuccess == 1 && timePeriod > 1800 ) %TBMAC2b
        dead = simulateEvent(TBparams.CatIVdeath*cat4deathScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIVageBrac, stateMat(catIVTBppl,Page), stateMat(catIVTBppl,Psex));
    else
        dead = simulateEvent(TBparams.CatIVdeath, 0, TBparams.TBCatIVageBrac, stateMat(catIVTBppl,Page), stateMat(catIVTBppl,Psex));
    end
    
    %%TB MAC ONLY: DEATH ADJUSTOR.  Half of these slated for death do not die, actually self cure.
    actuallyNotdead = rand(length(catIVTBppl),1) < 0.43;
    selfCures = (dead == 1 & actuallyNotdead == 1);
    dead = (dead == 1 & actuallyNotdead == 0);
    stateMat(catIVTBppl(selfCures),Phealth) = 0;  %self cures!
    stateMat(catIVTBppl(selfCures),Ptreatment) = 0;  %self cures and leaves trt!
    %end the TB mac death adjustor  

    %     nonSmokingDead = (stateMat(catIVTBppl(dead),Psmoking) == 0);
    %     deathsCounterMat(timePeriod,9) = sum(nonSmokingDead);
    deathsCounterMat(timePeriod,10) = sum(dead) - deathsCounterMat(timePeriod,9);  %smoking dead
    deathsCounterMat(timePeriod,17) = size( subsetterFunc(stateMat,catIVTBppl(dead),Phealth,[4]) ,1); %count their death
    stateMat(catIVTBppl(dead),Phealth) = 6;
    deathAgesTB = [deathAgesTB ;stateMat(catIVTBppl(dead),Page)];
    if timePeriod >= burnInTime
        runningDead(catIVTBppl(dead)) = 1;
        runningMDRDead(catIVTBppl(dead)) = 1;
    end
    
    %also, anyone who is 99 years old automatically dies, regardless of treatment status/health
    nintyNine = subsetterFuncInequality(stateMat, [-1],Page,99,Inf);
    stateMat(nintyNine,Phealth) = 6;
    deathAges = [deathAges ;stateMat(nintyNine,Page)];
    %    nonSmokingDead = (stateMat(nintyNine,Psmoking) == 0); %NOT DOING SMOKING ANYMORE
    %     %%COMMENTED OUT FOR SPEED
    %     deathsCounterMat(timePeriod,11) = sum(nonSmokingDead);  %smoking dead
    deathsCounterMat(timePeriod,12) = size(nintyNine,1) - deathsCounterMat(timePeriod,11);  %smoking dead
    deathsCounterMat(timePeriod, 19) = deathsCounterMat(timePeriod, 19) + size(intersect(nintyNine,detectedPpl_private),1);%count their death
    if timePeriod >= burnInTime
        runningDead(nintyNine) = 1;
    end
    
    %active sensTB related deaths: all untreated, catI, II, IV deaths - MDR deaths
    deathsCounterMat(timePeriod,18) = sum(deathsCounterMat(timePeriod,3:10)) - sum(deathsCounterMat(timePeriod,13:17));
    
    %for death rates debug
    % Count the number of people alive
    if durationYrs == 230
        if timePeriod >= burnInTime
            for fish = 1:numCurrentPpl
                if cohortFlagMat(fish,1) == 1
                    if stateMat(fish, Phealth) == 3 || stateMat(fish, Phealth) == 4  %count them in the actTBmatrix
                        aliveActTBAgesMat_after(timePeriod-burnInTime+1, MonSinceAct(fish), stateMat(fish,Psmoking)+1, ageFlagMat(fish)+1) ...
                            = aliveActTBAgesMat_after(timePeriod-burnInTime+1, MonSinceAct(fish), stateMat(fish,Psmoking)+1, ageFlagMat(fish)+1) + 1;
                    end
                end
            end
        end
    end
    % end for natural death rates debug
    
    if timePeriod >= burnInTime
        deathAges = [deathAges ;deathAgesTB];  %add the TB death ages list to the total death ages list;
        deathAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(deathAges,[0;10;15;20;30;40;50;60;70;80;90;110])';
        deathTBAgesMat_TBmac(timePeriod-burnInTime+1,:) = histc(deathAgesTB,[0;10;15;20;30;40;50;60;70;80;90;110])';
    end
    
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'TB deaths';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%%%%Updating flag%%%%%%%%%%%%%%%%%%
    updating = stateMat(1:numCurrentPpl,Phealth) < 5;  %5 is not yet born and 6 is dead
    if timePeriod >= burnInTime
        weightedCount = weightedSumMaker(TBparams.livingCost, 0, TBparams.livingCostAgeBrac, stateMat(updating,Page), stateMat(updating,Psex));
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (weightedCount); %incur cost of living if not dead
    end
    if ( LEbuilder == 1 && all(stateMat(:,Phealth) == 6) ) %everyone is dead
        timePeriod = totPeriods;  %skip to next cohort
    end
    %%%%%%%%%%%disease progression%%%%%%%%%%%%%%%%%%%%%
    
    % latent to active
    latentPeople = find(stateMat(1:numCurrentPpl,Phealth) == 1 | stateMat(1:numCurrentPpl,Phealth) == 2); % latent sens or MDR people
    
    %people recently infected
    latentLessT2Yrs = subsetterFuncInequality(stateMat,latentPeople,PtimePerInfect,timePeriod-TBparams.fastSlowActivationThreshold-1,timePeriod);
    
    %latentNonSmokers = subsetterFunc(stateMat,latentLessT2Yrs,Psmoking,[0]);  %nonsmokers only
    latentNonSmokers = latentLessT2Yrs;  %NOT DOING SMOKING ANYMORE
    ageAtInfection = max(0, ( double(stateMat(latentNonSmokers,Page)) - floor(double((timePeriod - stateMat(latentNonSmokers,PtimePerInfect))) /12) ) );
    incrOneMore =  (  mod((timePeriod - stateMat(latentNonSmokers,PtimePerInfect)),12) > stateMat(latentNonSmokers,PageIncr)  );
    ageAtInfection(incrOneMore) = max(0,ageAtInfection(incrOneMore) - 1);
    %OLD  activateTB = simulateEvent(TBparams.latentToActLessT2Yr, 0, TBparams.latentToActAgeBrac, ageAtInfection,[]);
    [trash,binned_ageAtInfection] = histc(ageAtInfection, TBparams.latentToActAgeBrac(:,1));
    [trash,binned_PtimePerInfect] = histc((timePeriod - stateMat(latentNonSmokers,PtimePerInfect)),TBparams.timeSinceInf);
    matIndex = binned_ageAtInfection + 4*binned_PtimePerInfect; %turn these into matrix indicies (remember index goes down the columns)
    probActivate = TBparams.latentToActLessT2Yr(matIndex);
    activateTB = rand(size(latentNonSmokers,1),1) <= probActivate;
    
    MdrCasesActivated = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2);  %count nonSmoker MDR activations
    MdrCasesEvolvedLatNowAct = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),PmdrEvolved) == 1); %count nonSmoker MDR activations who evolved
    %% COMMENTED OUT FOR SPEED
    activaCounterMat(timePeriod, 3) = sum(activateTB);  %count all the nonSmoker fast Activations
    activaCounterMat(timePeriod, 8) = sum(stateMat(latentNonSmokers(activateTB),Page) >= 15);
    %     activaCountSex(timePeriod, 1) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 1 & stateMat(latentNonSmokers(activateTB),Psex) == 1);  %count male nonSmoker DS fast activations
    %     activaCountSex(timePeriod, 9) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),Psex) == 1);  %count male nonSmoker MDR fast activations
    %     activaCountSex(timePeriod, 5) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 1 & stateMat(latentNonSmokers(activateTB),Psex) == 2);  %count female nonSmoker DS fast activations
    %     activaCountSex(timePeriod, 13) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),Psex) == 2);  %count female nonSmoker MDR fast activations
    %     latentCount(timePeriod, 1) = sum(stateMat(latentNonSmokers,Phealth) == 1 & stateMat(latentNonSmokers,Psex) == 1);  %count male nonSmoker DS fast latents
    %     latentCount(timePeriod, 9) = sum(stateMat(latentNonSmokers,Phealth) == 2 & stateMat(latentNonSmokers,Psex) == 1);  %count male nonSmoker MDR fast latents
    %     latentCount(timePeriod, 5) = sum(stateMat(latentNonSmokers,Phealth) == 1 & stateMat(latentNonSmokers,Psex) == 2);  %count female nonSmoker DS fast latents
    %     latentCount(timePeriod, 13) = sum(stateMat(latentNonSmokers,Phealth) == 2 & stateMat(latentNonSmokers,Psex) == 2);  %count female nonSmoker MDR fast latents
    stateMat(latentNonSmokers(activateTB), Phealth) = stateMat(latentNonSmokers(activateTB), Phealth) + 2;  %adding two just means latent sens people now have active sens, and latent MDR have active MDR
    stateMat(latentNonSmokers(activateTB), PactTime) = timePeriod;  %write time they activated
    fastActivators = latentNonSmokers(activateTB);  %DEBUGGING
    %NOT DOING SMOKING ANYMORE
    %     latentSmokers = subsetterFunc(stateMat,latentLessT2Yrs,Psmoking,[1]);  %smokers only
    %     activateTB = simulateEvent(TBparams.latentToActLessT2YrSmoker, 0, TBparams.latentToActAgeBrac, stateMat(latentSmokers,Page),[]);
    %     MdrCasesActivated = MdrCasesActivated + sum(stateMat(latentSmokers(activateTB),Phealth) == 2);  %count Smoker MDR activations
    %     MdrCasesEvolvedLatNowAct = MdrCasesEvolvedLatNowAct + sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),PmdrEvolved) == 1); %count Smoker MDR activations who evolved
    % %% COMMENTED OUT FOR SPEED
    %      activaCounterMat(timePeriod, 1) = size(latentLessT2Yrs,1); %count all the  less than 2 years latentTB
    %      activaCounterMat(timePeriod, 4) = sum(activateTB);  %count all the  smoker fast Activations
    % %     activaCountSex(timePeriod, 3) = sum(stateMat(latentSmokers(activateTB),Phealth) == 1 & stateMat(latentSmokers(activateTB),Psex) == 1);  %count male smoker DS fast activations
    % %     activaCountSex(timePeriod, 11) = sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),Psex) == 1);  %count male smoker MDR fast activations
    % %     activaCountSex(timePeriod, 7) = sum(stateMat(latentSmokers(activateTB),Phealth) == 1 & stateMat(latentSmokers(activateTB),Psex) == 2);  %count female smoker DS fast activations
    % %     activaCountSex(timePeriod, 15) = sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),Psex) == 2);  %count female smoker MDR fast activations
    % %     latentCount(timePeriod, 3) = sum(stateMat(latentSmokers,Phealth) == 1 & stateMat(latentSmokers,Psex) == 1);  %count male smoker DS fast latents
    % %     latentCount(timePeriod, 11) = sum(stateMat(latentSmokers,Phealth) == 2 & stateMat(latentSmokers,Psex) == 1);  %count male smoker MDR fast latents
    % %     latentCount(timePeriod, 7) = sum(stateMat(latentSmokers,Phealth) == 1 & stateMat(latentSmokers,Psex) == 2);  %count female smoker DS fast latents
    % %     latentCount(timePeriod, 15) = sum(stateMat(latentSmokers,Phealth) == 2 & stateMat(latentSmokers,Psex) == 2);  %count female smoker MDR fast latents
    %     stateMat(latentSmokers(activateTB), Phealth) = stateMat(latentSmokers(activateTB), Phealth) + 2;  %adding two just means latent sens people now have active sens, and latent MDR have active MDR
    %     stateMat(latentSmokers(activateTB), PactTime) = timePeriod;  %write time they activated
    %
    %END NOT DOING SMOKING ANYMORE
    %people infected more than 2 years ago
    latentGreaterT2Yrs = subsetterFuncInequality(stateMat,latentPeople,PtimePerInfect,1,timePeriod-TBparams.fastSlowActivationThreshold);
    
    %    latentNonSmokers =
    %    subsetterFunc(stateMat,latentGreaterT2Yrs,Psmoking,[0]);  %nonsmokers only
    latentNonSmokers = latentGreaterT2Yrs; % NOT DOING SMOKING ANYMORE
    ageAtInfection = max(0, ( double(stateMat(latentNonSmokers,Page)) - floor(double((timePeriod - stateMat(latentNonSmokers,PtimePerInfect))) /12) ) );
    incrOneMore =  (  mod((timePeriod - stateMat(latentNonSmokers,PtimePerInfect)),12) > stateMat(latentNonSmokers,PageIncr)  );
    ageAtInfection(incrOneMore) = max(0,ageAtInfection(incrOneMore) - 1);
    %OLD   activateTB = simulateEvent(TBparams.latentToActGreatT2Yr, 0, TBparams.latentToActAgeBrac, ageAtInfection,[]);
    
    timeSince2years = double( (timePeriod-24) - stateMat(latentNonSmokers,PtimePerInfect)) / 12;  %years from activtion, with year three = 1 since the equation indexes at 1.
    [trash,binned_slowPtimePerInfect] = histc(timeSince2years, TBparams.timeSinceAct_slowBins);
    probActivate = TBparams.latentToActGreatT2Yr(binned_slowPtimePerInfect)';
    activateTB = rand(size(latentNonSmokers,1),1) <= probActivate;
    
    MdrCasesActivated = MdrCasesActivated + sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2);  %count nonSmoker MDR activations
    MdrCasesEvolvedLatNowAct = MdrCasesEvolvedLatNowAct + sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),PmdrEvolved) == 1); %count nonSmoker MDR activations who evolved

    %% COMMENTED OUT FOR SPEED
    activaCounterMat(timePeriod, 5) = sum(activateTB);  %count all the  nonsmoker slow Activations
    %     activaCountSex(timePeriod, 2) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 1 & stateMat(latentNonSmokers(activateTB),Psex) == 1);  %count male nonSmoker DS slow activations
    %     activaCountSex(timePeriod, a10) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),Psex) == 1);  %count male nonSmoker MDR slow activations
    %     activaCountSex(timePeriod, 6) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 1 & stateMat(latentNonSmokers(activateTB),Psex) == 2);  %count female nonSmoker DS slow activations
    %     activaCountSex(timePeriod, 14) = sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),Psex) == 2);  %count female nonSmoker MDR slow activations
    %     latentCount(timePeriod, 2) = sum(stateMat(latentNonSmokers,Phealth) == 1 & stateMat(latentNonSmokers,Psex) == 1);  %count male nonSmoker DS slow latents
    %     latentCount(timePeriod, 10) = sum(stateMat(latentNonSmokers,Phealth) == 2 & stateMat(latentNonSmokers,Psex) == 1);  %count male nonSmoker MDR slow latents
    %     latentCount(timePeriod, 6) = sum(stateMat(latentNonSmokers,Phealth) == 1 & stateMat(latentNonSmokers,Psex) == 2);  %count female nonSmoker DS slow latents
    %     latentCount(timePeriod, 14) = sum(stateMat(latentNonSmokers,Phealth) == 2 & stateMat(latentNonSmokers,Psex) == 2);  %count female nonSmoker MDR slow latents
    stateMat(latentNonSmokers(activateTB), Phealth) = stateMat(latentNonSmokers(activateTB), Phealth) + 2;  %adding two just means latent sens people now have active sens, and latent MDR have active MDR
    stateMat(latentNonSmokers(activateTB), PactTime) = timePeriod;  %write time they activated
    slowActivators = latentNonSmokers(activateTB);  %DEBUGGING
    
    % old people activate
    latentGreaterT14Yrs = subsetterFuncInequality(stateMat,latentPeople,PtimePerInfect,1,timePeriod-(12*14));
    latentGreaterT14YrsOLD = subsetterFuncInequality(stateMat,latentGreaterT14Yrs,Page,30,101);
    cappedAge = double(min(stateMat(latentGreaterT14YrsOLD,Page),60));
    
    % activateTB = rand(size(latentGreaterT14YrsOLD,1),1) <= oldActSlope_cal + (  (1+oldActIntercept_cal).^((cappedAge-30)/5) - 1 )  ;
    activateTB = rand(size(latentGreaterT14YrsOLD,1),1) <= (oldActSlope_cal + (  (1+oldActIntercept_cal).^((cappedAge-30)/5) - 1 ))  ;
    oldActIntercept_cal*((1+oldActSlope_cal).^(stateMat(latentGreaterT14YrsOLD,Page)-30));
    ((( oldActSlope_cal )*(stateMat(latentGreaterT14YrsOLD,Page)-30)) - oldActIntercept_cal)  ;
    
    %write some stuff down
    MdrCasesActivated = MdrCasesActivated + sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2);  %count nonSmoker MDR activations
    MdrCasesEvolvedLatNowAct = MdrCasesEvolvedLatNowAct + sum(stateMat(latentNonSmokers(activateTB),Phealth) == 2 & stateMat(latentNonSmokers(activateTB),PmdrEvolved) == 1); %count nonSmoker MDR activations who evolved
    activaCounterMat(timePeriod, 5) = activaCounterMat(timePeriod, 5)+sum(activateTB);
    
    stateMat(latentGreaterT14YrsOLD(activateTB), Phealth) = stateMat(latentGreaterT14YrsOLD(activateTB), Phealth) + 2;  %adding two just means latent sens people now have active sens, and latent MDR have active MDR
    stateMat(latentGreaterT14YrsOLD(activateTB), PactTime) = timePeriod;  %write time they activated
    OLDslowActivators = latentGreaterT14YrsOLD(activateTB);  %DEBUGGING
    latentGreaterT14Yrs_prob =  oldActSlope_cal + (  (1+oldActIntercept_cal).^(([30:1:50]-30)/5) - 1 ) ; %DEBUGGING
    
    
    % NOT DOING SMOKING ANYMORE
    %     latentSmokers = subsetterFunc(stateMat,latentGreaterT2Yrs,Psmoking,[1]);  %smokers only
    %     activateTB = simulateEvent(TBparams.latentToActGreatT2YrSmoker, 0, TBparams.latentToActAgeBrac, stateMat(latentSmokers,Page),[]);
    %     MdrCasesActivated = MdrCasesActivated + sum(stateMat(latentSmokers(activateTB),Phealth) == 2);  %count Smoker MDR activations
    %     MdrCasesEvolvedLatNowAct = MdrCasesEvolvedLatNowAct + sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),PmdrEvolved) == 1); %count Smoker MDR activations who evolved
    %     activaCounterMat(timePeriod, 2) = size(latentGreaterT2Yrs,1); %count all the  more than 2 years latentTB
    %     activaCounterMat(timePeriod, 6) = sum(activateTB);  %count all the  smoker slow Activations
    % %% COMMENTED OUT FOR SPEED
    % %     activaCountSex(timePeriod, 4) = sum(stateMat(latentSmokers(activateTB),Phealth) == 1 & stateMat(latentSmokers(activateTB),Psex) == 1);  %count male smoker DS slow activations
    % %     activaCountSex(timePeriod, 12) = sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),Psex) == 1);  %count male smoker MDR slow activations
    % %     activaCountSex(timePeriod, 8) = sum(stateMat(latentSmokers(activateTB),Phealth) == 1 & stateMat(latentSmokers(activateTB),Psex) == 2);  %count female smoker DS slow activations
    % %     activaCountSex(timePeriod, 16) = sum(stateMat(latentSmokers(activateTB),Phealth) == 2 & stateMat(latentSmokers(activateTB),Psex) == 2);  %count female smoker MDR slow activations
    % %     latentCount(timePeriod, 4) = sum(stateMat(latentSmokers,Phealth) == 1 & stateMat(latentSmokers,Psex) == 1);  %count male smoker DS slow latents
    % %     latentCount(timePeriod, 12) = sum(stateMat(latentSmokers,Phealth) == 2 & stateMat(latentSmokers,Psex) == 1);  %count male smoker MDR slow latents
    % %     latentCount(timePeriod, 8) = sum(stateMat(latentSmokers,Phealth) == 1 & stateMat(latentSmokers,Psex) == 2);  %count female smoker DS slow latents
    % %     latentCount(timePeriod, 16) = sum(stateMat(latentSmokers,Phealth) == 2 & stateMat(latentSmokers,Psex) == 2);  %count female smoker MDR slow latents
    %     stateMat(latentSmokers(activateTB), Phealth) = stateMat(latentSmokers(activateTB), Phealth) + 2;  %adding two just means latent sens people now have active sens, and latent MDR have active MDR
    %     stateMat(latentSmokers(activateTB), PactTime) = timePeriod;  %write time they activated
    %END NOT DOING SMOKING ANYMORE
    
    %record and get willing to seek care
    activatedToday = subsetterFunc(stateMat,[-1],PactTime,timePeriod);
    activaCounterMat(timePeriod,7) = size(subsetterFuncInequality(stateMat,activatedToday,Page,0,14),1);  %count just activated kids
    
    %some of these people are willing to seek RNTCP care
    actToday_needAssignment = subsetterFunc(stateMat,activatedToday,PassignedSeekTrtmt,0);
    logical_seekRNTCP = simulateEvent(TBparams.probSeekRNTCPtreatment, 0, [0 Inf], stateMat(actToday_needAssignment,Page),[]);
    stateMat(actToday_needAssignment(logical_seekRNTCP),PseekTrtmt) = 1;
    %Commented out since we are going to do person-episodic RNTCP seeking
    %assignment (instead of by person)
    %    stateMat(activatedToday,PassignedSeekTrtmt) = 1;  %mark that they were assigned treatment-seeking status and never need to do so again
    %%COMMENTED OUT FOR SPEED
    %     activaToday(timePeriod, 1) = size(subsetterFunc(stateMat,actToday_needAssignment,Phealth,3),2);
    %     activaToday(timePeriod, 2) = size(subsetterFunc(stateMat,actToday_needAssignment,Phealth,4),2);
    
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'latent to active';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%% Treatment
    RNTCPrampUpPeriod = 10*12;
    if timePeriod >= burnInTime
        DotsPlusStartTime = burnInTime + RNTCPrampUpPeriod + 12; %starts one year after RNTCP full coverage;
        dotsPlusCoverage = 0;
        if (timePeriod > DotsPlusStartTime) %one year after rntcp full coverage
            if (timePeriod < DotsPlusStartTime + TBparams.DotsPlusRampUpPeriod )
                dotsPlusCoverage = TBparams.DotsPlusRampUpSequence(timePeriod-DotsPlusStartTime,1);
            end
        end
        
        % CAT I
        numCATI = size(find(stateMat(1:numCurrentPpl,Ptreatment) == 1),1); %record num ppl in catI
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + numCATI*(TBparams.monthlyNonmedPatientCosts + TBparams.catIdrugCost + TBparams.catIclinicCost); %incur costs from being in catI
        subsetCATI = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 1);

        subsetCATImales = subsetterFunc(stateMat, subsetCATI, Psex, [1]);  %count all the males
        subsetCATIdemog(1) = size(subsetCATImales, 1);
        subsetCATIdemog(2) = size(subsetterFuncInequality(stateMat, subsetCATImales, Page, 0, 45),1);  %count all the young males
        subsetCATIfemales = subsetterFunc(stateMat, subsetCATI, Psex, [2]);  %count all the females
        subsetCATIdemog(3) = size(subsetCATIfemales,1);
        subsetCATIdemog(4) = size(subsetterFuncInequality(stateMat, subsetCATIfemales, Page, 0, 45),1);  %count all the young females
        
        stateMat(subsetCATI,PtrtmtCounter) = stateMat(subsetCATI,PtrtmtCounter) + 1; % increment treatment counter
        defaults = simulateEvent(TBparams.CatIdefault, 0, TBparams.TBCatIageBrac, stateMat(subsetCATI,Page), stateMat(subsetCATI,Psex));
        CatIDefaulters = subsetCATI(defaults);  %these are defaulters
       
        if TBMAC2b_trtSuccess == 1 %TBMAC2b
            if timePeriod > 1800
                defaults = simulateEvent(TBparams.CatIdefault*cat12initScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIageBrac, stateMat(subsetCATI,Page), stateMat(subsetCATI,Psex));
                CatIDefaulters = subsetCATI(defaults);  %these are defaulters
            end
        end
        if TBMAC2a_initDefReduction == 1 %TBMAC2a
            if timePeriod > 1800
                catIinitialmo = subsetterFunc(stateMat, subsetCATI, PtrtmtCounter, 1);
                CatIDefaulters = setdiff(CatIDefaulters,catIinitialmo );  %these people cannot be marked to default yet
                defaults_initMo = simulateEvent(TBparams.CatIdefault*TbMac2_initDefaultDec(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIageBrac, stateMat(catIinitialmo,Page), stateMat(catIinitialmo,Psex));
                CatIDefaulters = sort([CatIDefaulters; catIinitialmo(defaults_initMo)]);
            end
        end
        
        numCATIdefaults = size(CatIDefaulters,1);  %record num defaults
        
        % Defaulted
        stateMat(CatIDefaulters,Ptreatment)     = 0; % defaulters go to no treatment
        stateMat(CatIDefaulters,PpastTreatment) = 2; % Mark as a defaulter
        stateMat(CatIDefaulters,PtrtmtCounter)  = 0; % Reset treatment counter
        if timePeriod >= burnInTime
            catInonSuccessTypes(timePeriod-burnInTime+1,1) = size(CatIDefaulters,1);  %record num of defaults
        end
        % Defaulted and got MDR tb
        subsetCATItb = subsetterFunc(stateMat, CatIDefaulters, Phealth, [3]);  %only if actually have TB (MDR ppl already have MDR)
        mdr = rand(length(subsetCATItb),1) <= TBparams.prbMDRafterDefault; % Got MDR
        stateMat(subsetCATItb(mdr),Phealth) = 4; % Set their health state to MDR TB
        MdrCasesEvolvedAct = sum(mdr);             %count the Mdr cases that evolved
        
        subsetCATI = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 1);  %now get non-defaulters
        
        % By 4 months DS latent TB people are cured
        fourMonthers = subsetterFunc(stateMat, subsetCATI, PtrtmtCounter, [4]);
        nonMDRfourMonthers_noTBnoMDR = subsetterFunc(stateMat, fourMonthers, Phealth, [1]); %latent DS TB ppl only
        stateMat(nonMDRfourMonthers_noTBnoMDR,Phealth) = 0; %treatment cures latent DS
        
        % Not defaulted and 6 months
        sixMonthers   = subsetterFunc(stateMat,subsetCATI,PtrtmtCounter,[6]);
        numCATIgraduating = size(sixMonthers,1); %record num ppl done with catI
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (numCATIgraduating*TBparams.SStestCost*2); % SS cost for all sixmonthers, 2 SS tests to establish cure
        
        activeDSppl = subsetterFunc(stateMat, sixMonthers, Phealth, [3]); %MDR people can't get cured
        curedLogi     = simulateEvent(TBparams.CatIsuccess, 0, TBparams.TBCatIageBrac, stateMat(activeDSppl,Page), stateMat(activeDSppl,Psex));
        curedNoMdrppl = activeDSppl(curedLogi);
        numCATIcured = size(curedNoMdrppl,1); %record num ppl 'cured'
        
        % Not defaulted and 6 months - "cured". May actually have latent TB.
        %some are actually cured.  Set all cured to healthy, then we will reset some to latent;
        stateMat(curedNoMdrppl,Phealth)        = 0;
        stateMat(curedNoMdrppl,Ptreatment)     = 0;
        stateMat(curedNoMdrppl,PtrtmtCounter)  = 0;
        stateMat(curedNoMdrppl,PpastTreatment) = 3;
        %some actually have latent TB still
        curedToLatent = rand(length(curedNoMdrppl),1) <= TBparams.prbCatIrelapse; %these people are actually not cured and have latent TB still
        stillLatent = curedNoMdrppl(curedToLatent);
        stateMat(stillLatent,Phealth) = 1;   %go to latent TB
        stateMat(stillLatent,PtimePerInfect) = timePeriod; %reset infection time, since newly enter latent pool.
        if timePeriod >= burnInTime
            catInonSuccessTypes(timePeriod-burnInTime+1, 2) = sum(curedToLatent);  %record latent "successes"
        end
        %some of those cured people actually get latent MDR
        curedToMDRLat = rand(length(curedNoMdrppl),1) <= TBparams.prbMDRtrtmtB4;  %some of these people get MDR latent;
        nowMDRlat = curedNoMdrppl(curedToMDRLat);
        stateMat(nowMDRlat,Phealth) = 2; %go to MDR latent
        stateMat(nowMDRlat,PtimePerInfect) = timePeriod; %reset infection time, since newly enter latent pool.
        MdrCasesEvolvedLat = size(nowMDRlat,1);  %count the Mdr cases that evolved
        stateMat(nowMDRlat, PmdrEvolved) = 1; %mark that mdr was evolved
        
        %some MDR ppl test neg for TB becuase ss test is not that good
        sixMonthersMDR = subsetterFunc(stateMat, sixMonthers, Phealth, [4]);
        testNeg = rand(length(sixMonthersMDR),1) <= (1-TBparams.SSsensit); %correct "cure" rate for MDR ppl is getting an incorrect sputum smear test (test say neg for TB)
        stateMat(sixMonthersMDR(testNeg),Ptreatment) = 0;
        stateMat(sixMonthersMDR(testNeg),PtrtmtCounter) = 0;
        stateMat(sixMonthersMDR(testNeg),PpastTreatment ) = 3;
        numCATIcured = numCATIcured + size(sixMonthersMDR(testNeg),1); %record num ppl 'cured'
        
        %for scenario where DOTS can cure MDR
        if TBparams.catIcuresMDR > 0
            curedMDR = rand(length(sixMonthersMDR),1) <= TBparams.catIcuresMDR;
            stateMat(sixMonthersMDR(curedMDR),Phealth) = 0;
        end
        
        %for non-active TB patients.  After this there will still be some non-activeTB patients in RNTCP if spec is not 1
        sixMonthers_noTB = subsetterFunc(stateMat, sixMonthers, Phealth, [0,1,2]);
        curedLogi_noTB = rand(length(sixMonthers_noTB),1) <= TBparams.SSspec;
        cured_noTB     = sixMonthers_noTB(curedLogi_noTB);
        stateMat(cured_noTB,Ptreatment)     = 0;  %released from treatment with no change in health status
        stateMat(cured_noTB,PtrtmtCounter)  = 0;
        stateMat(cured_noTB,PpastTreatment) = 3;
        
        % Not defaulted and 6 months - failed.  Go to CATII treatment.  (this includes all MDR ppl, in next paragraph rewrote some so they test pos)
        failed = setdiff(sixMonthers,sort(unique([curedNoMdrppl;sixMonthersMDR(testNeg);cured_noTB]))); %get sixMonthers who were not 'cured';
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(failed,1) * TBparams.DSTcost*dotsPlusCoverage);  %incur DST cost for MDR testing and adjust for coverage
        stateMat(failed,Ptreatment)    = 2;
        stateMat(failed,PtrtmtCounter) = 0;
        if timePeriod >= burnInTime
            catInonSuccessTypes(timePeriod-burnInTime+1, 3) = size(failed,1);  %record failures
        end
        % If not already going to test positive for MDR, overwrite MDR testing (note we are actually testing more people than just this)
        failedNotAlreadyInMDRtesting = subsetterFunc(stateMat, failed, PmdrTestingStatus, [0, 1, 2, 3]);
        stateMat(failedNotAlreadyInMDRtesting,PmdrTesting) = timePeriod; %test for MDR if failed
        stateMat(failedNotAlreadyInMDRtesting,PmdrTestingStatus) = stateMat(failedNotAlreadyInMDRtesting, Phealth);
        mdrTestPpl(timePeriod, 16) = size(subsetterFunc(stateMat, failedNotAlreadyInMDRtesting, Phealth, [0,1,2]),1);
        mdrTestPpl(timePeriod, 17) = size(subsetterFunc(stateMat, failedNotAlreadyInMDRtesting, Phealth, [3]),1);
        mdrTestPpl(timePeriod, 18) = size(subsetterFunc(stateMat, failedNotAlreadyInMDRtesting, Phealth, [4]),1);
        
        
        %=======================================%
        
        % CAT II
        numCATII = size(find(stateMat(1:numCurrentPpl,Ptreatment) == 2),1); %record num ppl in catII
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + numCATII*(TBparams.monthlyNonmedPatientCosts + TBparams.catIIdrugCost + TBparams.catIIclinicCost );  %incur costs from being in catII
        subsetCATII = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 2);
        if timePeriod >= burnInTime
            subsetCATIImales = subsetterFunc(stateMat, subsetCATII, Psex, [1]);  %count all the males
            subsetCATIIdemog(1) = size(subsetCATIImales,1);
            subsetCATIIdemog(2) = size(subsetterFuncInequality(stateMat, subsetCATIImales, Page, 0, 45),1);  %count all the young males
            subsetCATIIfemales = subsetterFunc(stateMat, subsetCATII, Psex, [2]);  %count all the females
            subsetCATIIdemog(3) = size(subsetCATIIfemales,1);
            subsetCATIIdemog(4) = size(subsetterFuncInequality(stateMat, subsetCATIIfemales, Page, 0, 45),1);  %count all the young females
        end
        
        stateMat(subsetCATII,PtrtmtCounter) = stateMat(subsetCATII,PtrtmtCounter) + 1; %increment treatment counter
        defaults = simulateEvent(TBparams.CatIIdefault, 0, TBparams.TBCatIIageBrac, stateMat(subsetCATII,Page), stateMat(subsetCATII,Psex));
        
        CatIIDefaulters = subsetCATII(defaults);  %these are defaulters

        if TBMAC2b_trtSuccess == 1
            if timePeriod > 1800
                defaults = simulateEvent(TBparams.CatIIdefault*cat12initScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIIageBrac, stateMat(subsetCATII,Page), stateMat(subsetCATII,Psex));
                CatIIDefaulters = subsetCATII(defaults);
            end
        end
        if TBMAC2a_initDefReduction == 1
            if timePeriod > 1800
                catIIinitialmo = subsetterFunc(stateMat, subsetCATII, PtrtmtCounter, 1);
                CatIIDefaulters = setdiff(CatIIDefaulters,catIIinitialmo );  %these people cannot be marked to default yet
                defaults_initMo = simulateEvent(TBparams.CatIIdefault*TbMac2_initDefaultDec(min(1884-1800,timePeriod-1800))*TBparams.TBmacMultiplier, 0, TBparams.TBCatIIageBrac, stateMat(catIIinitialmo,Page), stateMat(catIIinitialmo,Psex));
                CatIIDefaulters = sort([CatIIDefaulters; catIIinitialmo(defaults_initMo)]);
            end
        end
        
        numCATIIdefaults = sum(defaults);  %record num defaults
        
        % Defaulted
        stateMat(CatIIDefaulters,Ptreatment)     = 0; % defaulters go to no treatment
        stateMat(CatIIDefaulters,PpastTreatment) = 2;
        stateMat(CatIIDefaulters,PtrtmtCounter)  = 0;
        % Defaulted and got MDR tb
        DScatIIDefaulters = subsetterFunc(stateMat, CatIIDefaulters, Phealth, [3]);
        mdr = rand(length(DScatIIDefaulters),1) <= TBparams.prbMDRafterDefault;
        stateMat(DScatIIDefaulters(mdr), Phealth) = 4;        % got MDR tb
        MdrCasesEvolvedAct = MdrCasesEvolvedAct + sum(mdr); % count the Mdr cases that evolved
        
        subsetCATII = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 2);  %now only want nondefaulters
        
        % By 4 months DS latent TB people are cured
        fourMonthers = subsetterFunc(stateMat, subsetCATII, PtrtmtCounter, [4]);
        nonMDRfourMonthers_noTBnoMDR = subsetterFunc(stateMat, fourMonthers, Phealth, [1]); %latent DS TB ppl only
        stateMat(nonMDRfourMonthers_noTBnoMDR,Phealth) = 0; %treatment cures latent DS
        
        % 4 months: all smear pos get tested for MDR
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(fourMonthers,1)*TBparams.SStestCost); % SS cost for all fourMonthers
        nonMDRfourMonthers = subsetterFunc(stateMat, fourMonthers, Phealth, [3]); %sensTB ppl only
        testPosAt4Mon      = rand(length(nonMDRfourMonthers),1) <= TBparams.testPos4MonCatII;
        % all MDR people would be smear pos, but SS test isn't perfect
        MDRfourMonthers = subsetterFunc(stateMat, fourMonthers, Phealth, [4]); %mdrTB ppl only
        testPosMDR      = rand(length(MDRfourMonthers),1) <= TBparams.SSsensit;
        % all noTB people would be smear neg, but test may have imperfect spec
        noTBfourMonthers = subsetterFunc(stateMat, fourMonthers, Phealth, [0,1,2]); %mdrTB ppl only
        testPosNoTB      = rand(length(noTBfourMonthers),1) <= 1-TBparams.SSspec;
        % If not already going to test positive for MDR, overwrite MDR testing (note we are actually testing more people than just this)
        allMDRtestEligible     = unique(sort([nonMDRfourMonthers(testPosAt4Mon); MDRfourMonthers(testPosMDR);noTBfourMonthers(testPosNoTB)]));
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(allMDRtestEligible,1) * TBparams.DSTcost*dotsPlusCoverage);  %incur DST cost for MDR testing and adjust for coverage
        notAlreadyInMDRtesting = subsetterFunc(stateMat, allMDRtestEligible, PmdrTestingStatus, [0, 1, 2, 3]);
        stateMat(notAlreadyInMDRtesting,PmdrTesting)       = timePeriod;
        stateMat(notAlreadyInMDRtesting,PmdrTestingStatus) = stateMat(notAlreadyInMDRtesting, Phealth);
        mdrTestPpl(timePeriod, 19) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [0,1,2]),1);
        mdrTestPpl(timePeriod, 20) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [3]),1);
        mdrTestPpl(timePeriod, 21) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [4]),1);
        
        % 8 months
        eightMonthers = subsetterFunc(stateMat,subsetCATII,PtrtmtCounter,[8]); %month when treatment ends
        numCATIIgraduating = size(eightMonthers,1); %record num ppl done with treatment
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 1) = size(eightMonthers,1);  %record num8monthers
        end
        
        curedNoMdrppl = subsetterFunc(stateMat, eightMonthers, Phealth, [3]);   %no cure if mdr.  Only cured if sensTB
        curedLogi     = simulateEvent(TBparams.CatIIsuccess, 0, TBparams.TBCatIIageBrac, stateMat(curedNoMdrppl,Page), stateMat(curedNoMdrppl,Psex));
        curedNoMdrppl = curedNoMdrppl(curedLogi);
        numCATIIcured = size(curedNoMdrppl,1); %record num ppl 'cured'
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 2) = size(curedNoMdrppl,1);  %record DS 8monthers
        end
        
        % 8 months - cured. Go to healthy
        stateMat(curedNoMdrppl,Phealth)        = 0;
        stateMat(curedNoMdrppl,Ptreatment)     = 0;
        stateMat(curedNoMdrppl,PtrtmtCounter)  = 0;
        stateMat(curedNoMdrppl,PpastTreatment) = 3;
        
        %some actually have latent TB still
        curedToLatent = rand(length(curedNoMdrppl),1) <= TBparams.prbCatIIrelapse; %these people are actually not cured and have latent TB still
        stillLatent   = curedNoMdrppl(curedToLatent);
        stateMat(stillLatent,Phealth)        = 1;          %go to latent TB
        stateMat(stillLatent,PtimePerInfect) = timePeriod; %reset infection time, since newly enter latent pool.
        %some of those cured people actually get latent MDR
        curedToMDRLat = rand(length(curedNoMdrppl),1) <= TBparams.prbMDRtrtmtB4;  %some of these people get MDR latent;
        nowMDRlat     = curedNoMdrppl(curedToMDRLat);
        stateMat(nowMDRlat, Phealth)        = 2;          %go to MDR latent
        stateMat(nowMDRlat, PtimePerInfect) = timePeriod; %reset infection time, since newly enter latent pool.
        MdrCasesEvolvedLat = MdrCasesEvolvedLat + size(nowMDRlat,1);  %count the Mdr cases that evolved
        stateMat(nowMDRlat, PmdrEvolved) = 1; %mark that mdr was evolved
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 3) = size(nowMDRlat,1);  %record lat mdr acquisition
        end
        
        %some MDR ppl test neg for TB becuase ss test is not that good (above code assumes all MDR ppl do not get 'cured')
        eightMonthersMDR = subsetterFunc(stateMat, eightMonthers, Phealth, [4]);
        testNeg = rand(length(eightMonthersMDR),1) <= (1-TBparams.SSsensit); % "cure" rate for MDR ppl is getting an incorrect sputum smear test (test say neg for TB)
        stateMat(eightMonthersMDR(testNeg),Ptreatment)     = 0;
        stateMat(eightMonthersMDR(testNeg),PtrtmtCounter)  = 0;
        stateMat(eightMonthersMDR(testNeg),PpastTreatment) = 1;
        numCATIIcured = numCATIIcured + size(eightMonthersMDR(testNeg),1);  %record num of these people as 'cured'
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 4) = size(eightMonthersMDR(testNeg),1);  %record MDR 'cured'
        end
        
        %for scenario where DOTS can cure MDR
        if TBparams.catIcuresMDR > 0
            curedMDR = rand(length(eightMonthersMDR),1) <= TBparams.catIIcuresMDR;
            stateMat(eightMonthersMDR(curedMDR),Phealth) = 0;
        end
        
        %for non-active TB patients.  After this there will still be some non-activeTB patients in RNTCP if spec is not 1
        eightMonthers_noTB = subsetterFunc(stateMat, eightMonthers, Phealth, [0,1,2]);
        curedLogi_noTB = rand(length(eightMonthers_noTB),1) <= TBparams.SSspec;
        cured_noTB     = eightMonthers_noTB(curedLogi_noTB);
        stateMat(cured_noTB,Ptreatment)     = 0;  %released from treatment with no change in health status
        stateMat(cured_noTB,PtrtmtCounter)  = 0;
        stateMat(cured_noTB,PpastTreatment) = 3;
        
        % 8 months - failed. Drop out of treatment.
        failed = setdiff(eightMonthers,sort(unique([curedNoMdrppl;eightMonthersMDR(testNeg);cured_noTB])));  % eightMonthers who were not successfully "cured" ;
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(failed,1) * TBparams.DSTcost*dotsPlusCoverage);  %incur DST cost for all failed
        stateMat(failed,Ptreatment)     = 0;
        stateMat(failed,PtrtmtCounter)  = 0;
        stateMat(failed,PpastTreatment) = 1;
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 5) = size(failed,1);  %record failed
        end
        % Some DS failed people get MDR
        DSfailed = subsetterFunc(stateMat, failed, Phealth, [3]);
        newMDR = rand(length(DSfailed),1) <= TBparams.prbMDRafterFailure;  %might get mdr tb
        stateMat(DSfailed(newMDR),Phealth) = 4;
        MdrCasesEvolvedAct = MdrCasesEvolvedAct + sum(newMDR);  %count the Mdr cases that evolved
        if timePeriod >= burnInTime
            catIInonSuccessTypes(timePeriod-burnInTime+1, 6) = sum(newMDR);  %record failed got MDR
        end
        % If not already going to test positive for MDR, overwrite MDR testing (note we are actually testing more people than just this)
        notAlreadyInMDRtesting = subsetterFunc(stateMat, failed, PmdrTestingStatus, [0,1,2,3]);
        stateMat(notAlreadyInMDRtesting,PmdrTesting) = timePeriod;
        stateMat(notAlreadyInMDRtesting,PmdrTestingStatus) = stateMat(notAlreadyInMDRtesting, Phealth);
        mdrTestPpl(timePeriod, 22) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [0,1,2]),1);
        mdrTestPpl(timePeriod, 23) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [3]),1);
        mdrTestPpl(timePeriod, 24) = size(subsetterFunc(stateMat, notAlreadyInMDRtesting, Phealth, [4]),1);
        
        %write some treatment demogs
        if timePeriod >= burnInTime
            treatmentDemog(timePeriod- burnInTime + 1, :) = [subsetCATIdemog(1),subsetCATIdemog(2), subsetCATIdemog(3), subsetCATIdemog(4), subsetCATIIdemog(1),subsetCATIIdemog(2), subsetCATIIdemog(3),subsetCATIIdemog(4) , numCurrentPpl ];
        end
        %=======================================%
        
        %CAT IV
        numCATIV = size(find(stateMat(1:numCurrentPpl,Ptreatment) == 4),1);  %record the num people in cat4
        subsetCATIV = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 4);
        
        stateMat(subsetCATIV,PtrtmtCounter) = stateMat(subsetCATIV,PtrtmtCounter) + 1; %increment treatment counter
        numCATIV_IP = size(find(stateMat(1:numCurrentPpl,Ptreatment) == 4 & stateMat(1:numCurrentPpl,PtrtmtCounter) <= 6),1);  %num of IP people in catIV
        numCATIV_CP = numCATIV - numCATIV_IP;
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + numCATIV_IP*(TBparams.monthlyNonmedPatientCosts + TBparams.catIVdrugIPcost + TBparams.catIVclinicCost);  %incur costs from being in catIV
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + numCATIV_CP*(TBparams.monthlyNonmedPatientCosts + TBparams.catIVdrugCPcost + TBparams.catIVclinicCost);  %incur costs from being in catIV
        
        defaults = simulateEvent(TBparams.CatIVdefault, 0, TBparams.TBCatIVageBrac, stateMat(subsetCATIV,Page), stateMat(subsetCATIV,Psex));
        CATIVDefaulters = subsetCATIV(defaults);  %these are defaulters
        
        if TBMAC2c_MDRtrtSuccess == 1
            if timePeriod > 1800
                defaults = simulateEvent(TBparams.CatIVdefault*cat4initScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIVageBrac, stateMat(subsetCATIV,Page), stateMat(subsetCATIV,Psex));
                CATIVDefaulters = subsetCATIV(defaults);  %these are defaulters
            end
        end
        if TBMAC2a_initDefReduction == 1
            if timePeriod > 1800
                catIVinitialmo = subsetterFunc(stateMat, subsetCATIV, PtrtmtCounter, 1);
                CATIVDefaulters = setdiff(CATIVDefaulters,catIVinitialmo );  %these people cannot be marked to default yet
                defaults_initMo = simulateEvent(TBparams.CatIVdefault*TbMac2_MDRinitDefaultDec(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIVageBrac, stateMat(catIVinitialmo,Page), stateMat(catIVinitialmo,Psex));
                CATIVDefaulters = sort([CATIVDefaulters; catIVinitialmo(defaults_initMo)]);
            end
        end
        
        %defaulted
        stateMat(CATIVDefaulters,Ptreatment)     = 0; % defaulters go to no treatment
        stateMat(CATIVDefaulters,PpastTreatment) = 2;
        stateMat(CATIVDefaulters,PtrtmtCounter)  = 0;
        
        % DS ppl Defaulted and got MDR tb
        DScatIVDefaulters = subsetterFunc(stateMat, CATIVDefaulters, Phealth, [3]);
        mdr = rand(length(DScatIVDefaulters),1) <= TBparams.prbMDRafterDefault;
        stateMat(DScatIVDefaulters(mdr), Phealth) = 4;        % got MDR tb
        MdrCasesEvolvedAct = MdrCasesEvolvedAct + sum(mdr); % count the Mdr cases that evolved
        
        subsetCATIV = find(updating & stateMat(1:numCurrentPpl,Ptreatment) == 4);  %re-get ppl still in treatment
        
        % By 4 months DS latent TB people are cured -- should we add latent MDR too? No, latent MDR gets cured at same rate and time as active MDR
        fourMonthers = subsetterFunc(stateMat, subsetCATIV, PtrtmtCounter, [4]);
        nonMDRfourMonthers_noTBnoMDR = subsetterFunc(stateMat, fourMonthers, Phealth, [0,1]); %no MDR non-active TB ppl only
        stateMat(nonMDRfourMonthers_noTBnoMDR,Phealth) = 0; %treatment cures latent DS
        
        % 24 months
        twenty4Monthers = subsetterFunc(stateMat, subsetCATIV, PtrtmtCounter, [24]); %month when treatment ends
        numCATIVgraduating = size(twenty4Monthers,1); %record num ppl done with catI
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(twenty4Monthers,1) * TBparams.DSTcost*dotsPlusCoverage);  %incur DST cost for all twenty2monthers
        twenty4Monthers_withMDR = subsetterFunc(stateMat, twenty4Monthers, Phealth, [2,4]); %MDR people might get cured
        cured = simulateEvent(TBparams.CatIVsuccess, 0, TBparams.TBCatIVageBrac, stateMat(twenty4Monthers_withMDR,Page), stateMat(twenty4Monthers_withMDR,Psex));
        if TBMAC2c_adv == 1
            if timePeriod > 1800
                cured = simulateEvent(TBparams.CatIVsuccess*cat4failureScalar(min(1884-1800,timePeriod-1800)), 0, TBparams.TBCatIVageBrac, stateMat(twenty4Monthers_withMDR,Page), stateMat(twenty4Monthers_withMDR,Psex));
            end
        end
        MDR_cured = twenty4Monthers_withMDR(cured);
        failed = twenty4Monthers_withMDR(~cured);
        
        % 24 months - cured. Go to healthy or latent.
        curedTwenty4Monthers = setdiff(twenty4Monthers,failed);  %everyone except MDR failures
        stateMat(curedTwenty4Monthers,Phealth)        = 0;
        stateMat(curedTwenty4Monthers,Ptreatment)     = 0;
        stateMat(curedTwenty4Monthers,PtrtmtCounter)  = 0;
        stateMat(curedTwenty4Monthers,PpastTreatment) = 3;
        stateMat(curedTwenty4Monthers,PmdrEvolved) = 0;
        
        %some actually have latent TB still
        curedToLatent = rand(length(MDR_cured),1) <= TBparams.prbCatIVrelapse; %these people are actually not cured and have latent TB still
        stillLatent   = MDR_cured(curedToLatent);
        stateMat(stillLatent, Phealth) = 2; %go to MDR latent
        
        % 24 months - failed. Released without further treatment
        stateMat(failed,Ptreatment)     = 0; %drop out of treatment
        stateMat(failed,PtrtmtCounter)  = 0;
        stateMat(failed,PpastTreatment) = 1;
        
        %record the DOTS numbers
        if timePeriod >= burnInTime
            DOTSperformanceMeas(timePeriod-burnInTime+1,:) = [numCATI, numCATIdefaults, numCATIgraduating, numCATIcured, numCATII, numCATIIdefaults, numCATIIgraduating, numCATIIcured, numCATIV, numCATIVgraduating, size(MDR_cured,1) ];
            actTBforCounts = subsetterFunc(stateMat, [-1], Phealth, [3,4]);
            pastTreatment(timePeriod-burnInTime+1, 1) = sum(stateMat(actTBforCounts, PpastTreatment) == 0);
            pastTreatment(timePeriod-burnInTime+1, 2) = sum(stateMat(actTBforCounts, PpastTreatment) == 1);
            pastTreatment(timePeriod-burnInTime+1, 3) = sum(stateMat(actTBforCounts, PpastTreatment) == 2);
            pastTreatment(timePeriod-burnInTime+1, 4) = sum(stateMat(actTBforCounts, PpastTreatment) == 3);
            pastTreatment(timePeriod-burnInTime+1, 5) = sum(stateMat(actTBforCounts, PprivTrtCount) ~= 0);  %has ever been on private trt in past
        end
        
    end %end "if timePeriod >= burnInTime"
    
    %=======================================%
    
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'treatment';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%start prepping for transmission stuff%%%%%%%%%%%
    if LEbuilder == 0
        ageBracTran  = ageBracketMaker(TBparams.modTransMatAgeBrac, stateMat(1:numCurrentPpl,Page));
        
        %prep transmission background vars
        %get sick people in each age bracket so you know transmission probs
        transmAgeBrac_sens = zeros(1,15);
        transmAgeBrac_mdr  = zeros(1,15);
        totPplInAgeBrac    = zeros(1,15);
        
        %summing by individual since it is faster than vector-element sum
        for i = 1:numCurrentPpl
            if (stateMat(i,Phealth) < 5)
                ageBracThisPerson =  ageBracTran(i) ; %if the person is in the age brac and not dead/unborn
                totPplInAgeBrac(1,ageBracThisPerson) = totPplInAgeBrac(1,ageBracThisPerson) + 1;
                if (stateMat(i,Phealth) == 3 && stateMat(i,Ptreatment) == 0)  %if has active sensTB & not being treated
                    transmAgeBrac_sens(1,ageBracThisPerson) = transmAgeBrac_sens(1,ageBracThisPerson) + 1;
                elseif (stateMat(i,Phealth) == 4 && stateMat(i,Ptreatment) ~= 4) %if has active mdrTB & not being treated for MDR
                    transmAgeBrac_mdr(1,ageBracThisPerson) = transmAgeBrac_mdr(1,ageBracThisPerson) + 1;
                end
            end
        end
        fracSensTb = transmAgeBrac_sens ./ totPplInAgeBrac;
        fracMdrTb  = transmAgeBrac_mdr  ./ totPplInAgeBrac;
        
        %for the minor revision scenario where MDR people in DOTS are less transmissible
        if TBparams.MDRinDOTS2lessTransmissible > 0
            MDRdotsPplInAgeBrac    = zeros(1,15);
            MDRpeople = subsetterFunc(stateMat,-1,Phealth,4);
            MDRinDOTS = subsetterFunc(stateMat,MDRpeople,Ptreatment,[2]);
            for i = 1:size(MDRinDOTS,1)
                ageBracThisPerson =  ageBracTran(i) ; %if the person is in the age brac and not dead/unborn
                MDRdotsPplInAgeBrac(1,ageBracThisPerson) = MDRdotsPplInAgeBrac(1,ageBracThisPerson) + 1;
            end
            fracMdrTb  = (transmAgeBrac_mdr-((1-TBparams.MDRinDOTS2lessTransmissible)*MDRdotsPplInAgeBrac)) ./ totPplInAgeBrac;
        end
        
        %for the scenario where kids don't transmit AT ALL
        fracSensTb(1:3) = TBparams.fracKidsSSpos*fracSensTb(1:3);
        fracMdrTb(1:3) = TBparams.fracKidsSSpos*fracMdrTb(1:3);
        
        % if dividing by zero, not want any transmission
        fracSensTb(totPplInAgeBrac == 0) = 0;
        fracMdrTb( totPplInAgeBrac == 0) = 0;
        
        %make rates (vector over age)
        probSensTB = TBparams.modMonRateMat * fracSensTb';
        probMdrTB  = TBparams.modMonRateMat * fracMdrTb';
        
        %%%%%%%%%%%%%%%%%%%%%%make transmission risks%%%%%%%%%%%%%%%%%%
        riskSensTB = 1-exp(-probSensTB(ageBracTran));
        riskMdrTB =  1-exp(-probMdrTB(ageBracTran));
        
        %timer
        timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
        timeLoopTimerName{timeLoopIndex} = 'prep transmission';
        timeLoopIndex = timeLoopIndex + 1;
        timerStart = cputime;
        
        %%%%%%%%%%%%%%%%Transmission%%%%%%%%%%%%%%
        susecptible = find(stateMat(1:numCurrentPpl,Phealth) == 0 );  %healthy, susceptible ppl
        susecptibleWithLat = find(stateMat(1:numCurrentPpl,Phealth) == 0 | stateMat(1:numCurrentPpl,Phealth) == 1 | stateMat(1:numCurrentPpl,Phealth) == 2 );  %for making effective contact rate
        %do sensTB before mdr -- sensTB has protective effect from getting Mdr
        
        %sensTB transmission
        SensTransmitted = rand(length(susecptible),1) <= riskSensTB(susecptible);
        hadTBbefore_ds = sum(stateMat(susecptible(SensTransmitted),PtimePerInfect) > 0);  %num ppl had TB before
        stateMat(susecptible(SensTransmitted),Phealth)        = 1;           %got latent sensTB
        stateMat(susecptible(SensTransmitted),PtimePerInfect) = timePeriod;  %mark when got infected
        DSCasesTransmitted = size(susecptible(SensTransmitted),1);           %count how many people got transmitted to
        
        %MDR transmission
        remainingSusecptible = find(stateMat(1:numCurrentPpl,Phealth) == 0 );
        MDRtransmitted       = rand(length(remainingSusecptible),1) <= TBparams.MdrRelFitness*riskMdrTB(remainingSusecptible);
        hadTBbefore_mdr = sum(stateMat(remainingSusecptible(MDRtransmitted),PtimePerInfect) > 0);  %num ppl had TB before
        stateMat(remainingSusecptible(MDRtransmitted), Phealth)       = 2;           %got MDR, latent
        stateMat(remainingSusecptible(MDRtransmitted),PtimePerInfect) = timePeriod;  %mark when got infected
        MdrCasesTransmitted = size(remainingSusecptible(MDRtransmitted),1);          %count how many people got transmitted to
        
        %%COMMENTED OUT FOR SPEED
        %     %count number of kids between 1 and 9 have infection
        %     kids = subsetterFunc(stateMat,[-1],Page,[1,9]); %kids
        %     parfor healthstateC = 0:4
        %         kids_TBstate = subsetterFunc(stateMat,kids,Phealth,[healthstateC]);
        %         kids_TBstateCount(timePeriod,healthstateC+1) = size(kids_TBstate,1);
        %     end
        %     %mean and median age of kids without active disease
        %     kids_noActTB = subsetterFunc(stateMat,kids,Phealth,[0,1,2]); %kids without active disease
        %     kids_TBstateCount(timePeriod,6) = median(single(stateMat(kids_noActTB, Page)));
        %     kids_TBstateCount(timePeriod,7) = mean(stateMat(kids_noActTB, Page));
        %
        %calculate force of infection, rather, the "effective contact rate"
        %     if mod(timePeriod,1)==0
        %         %forceOfInfection((timePeriod/12)/10,1) = forceOfInfectionCalculator(totPplInAgeBrac, transmAgeBrac_sens, transmAgeBrac_mdr);
        %         %forceOfInfection((timePeriod/12)/1,1) =  (DSCasesTransmitted + MdrCasesTransmitted)/ (sum(transmAgeBrac_sens) + sum(transmAgeBrac_mdr));
        %         forceOfInfection(timePeriod,1) = forceOfInfectionCalculator(susecptibleWithLat, riskSensTB, riskMdrTB, sum(transmAgeBrac_sens) + sum(transmAgeBrac_mdr)  );
        %     end
        
    end
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'transmission';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    
    %%%%%%%%%%%%%%%%NORMAL DIAGNOSIS%%%%%%%%%%%%%%%%%%%%%
    if (timePeriod >= burnInTime)
        %TBprevalence
        if TBparams.SSspec < 1
            activeAllPpl     = subsetterFunc(stateMat,-1,Ptreatment,[0]); % people not on treatment
        else
            activeTBppl_all = subsetterFunc(stateMat,-1,Phealth,[3,4]); %active people
            activeAllPpl     = subsetterFunc(stateMat,activeTBppl_all,Ptreatment,[0]); %active people not on treatment
        end
        
        %DOTS ramp up through coverage
        coverage = 0;
        if (timePeriod > burnInTime && timePeriod <= burnInTime + RNTCPrampUpPeriod)
            coverage = TBparams.DotsCoverageSequence(timePeriod-burnInTime,1);
        elseif timePeriod > (burnInTime+RNTCPrampUpPeriod)
            coverage = TBparams.DotsFullCoverage;
        end
        
        logical_coveredActTBppl = rand(length(activeAllPpl),1) <= coverage;
        coveredAllPpl_maySeek = activeAllPpl(logical_coveredActTBppl);
        
        %some people willing to seek RNTCP care
        coveredAllPpl = subsetterFunc(stateMat,coveredAllPpl_maySeek,PseekTrtmt,[1]);  %subset to those willing to seek care
        
        %identify TB knowledge for each person
        knowledgeAgeBrac = ageBracketMaker(TBparams.uptakeAgeBracs, stateMat(1:numCurrentPpl,Page));
        urbanKnowl = [TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,1), TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,2)];
        TBknowledge = -1*ones(numCurrentPpl,1);
        for ii = 1:size(coveredAllPpl,1)
            TBknowledge(coveredAllPpl(ii,1),1) = urbanKnowl(coveredAllPpl(ii,1),stateMat(coveredAllPpl(ii,1),Psex) );
        end
        
        %increase the prob for people with prior treatment.  Careful when
        %doing this that the ultimate prob of diagnosis does not exceed 1
        priorTreatBoost = ones(numCurrentPpl, 1);
        retreatmentTest = subsetterFunc(stateMat,coveredAllPpl,PpastTreatment,[1,2,3]);
        priorTreatBoost(retreatmentTest) = TBparams.priorTreatBoostFactor;  % max is 46.6483183
        catIIuptakeKnowledgeAges = TBparams.catIIuptakeKnowledge(knowledgeAgeBrac,:);
        for ii = 1:size(retreatmentTest,1)
            TBknowledge(retreatmentTest(ii)) = catIIuptakeKnowledgeAges(retreatmentTest(ii), stateMat(retreatmentTest(ii), Psex));  %prior treatment people must already know about treatment (ie, 1), averageKnowledgeOfTB is to compensate for 1/averageKnowledgeOfTB in the prTestedGivenTB param
        end
        
        if (timePeriod == burnInTime + 120) || (timePeriod == 1801)
            testingCoefficient
            TBparams.averageKnowledgeOfTB
            aveUptake_cal
            TBparams.prTestedGivenTB*TBparams.SSsensit * TBparams.uptakeUrbanKnowledge
            TBparams.prTestedGivenTB*(1-TBparams.SSspec)* TBparams.uptakeUrbanKnowledge
            TBparams.prTestedGivenTB*TBparams.SSsensit* TBparams.catIIuptakeKnowledge .* TBparams.priorTreatBoostFactor
            TBparams.prTestedGivenTB*(1-TBparams.SSspec)* TBparams.catIIuptakeKnowledge .* TBparams.priorTreatBoostFactor
        end
        
        %get the proper detection prob
        sensitUnique = ones(numCurrentPpl, 1)*(1-TBparams.SSspec);
        activeTBppl_all = subsetterFunc(stateMat,[-1],Phealth,[3,4]); %active people
        sensitUnique(activeTBppl_all) = TBparams.SSsensit;
        
        if TBMAC3b_XpertReplacesSmear_SSsensUp == 1
            if timePeriod > 1800
                activeCovered_tbmac = rand(size(activeTBppl_all,1),1) < GeneXcoverageMat(min(1848-1800, timePeriod-1800));
                sensitUnique(activeTBppl_all(activeCovered_tbmac)) = 0.65; %replace covered people with higher ssensit
            end
        end
        
        %diagnose people
        detectionProb = TBparams.prTestedGivenTB*sensitUnique.*TBknowledge.*priorTreatBoost;
        
        if TBMAC1a_increaseTrt == 1  %TB MAC scenario that increases uptake
            if (timePeriod > 1800)
                detectionProb = detectionProb * TBmac1ScalingFac(min(1884-1800, timePeriod-1800));  %jan 2016 to dec 2022                
            end
        end
        
        %cut off the ends on detection prob
        if any(detectionProb > 1)
            detectionProb(detectionProb > 1) = 1;
        end
        if any(detectionProb < 0)
            detectionProb(detectionProb < 0) = 0;
        end
        
        %DEBUGGING
        detectionProbDEBUGGING = TBparams.prTestedGivenTB*TBknowledge.*priorTreatBoost;
        canBeDetectedDEBUGGING = rand(numCurrentPpl,1) < detectionProbDEBUGGING;
        detectedDEBUGGING = rand(numCurrentPpl,1) < sensitUnique;
        indexVecForAllPpl = [1:numCurrentPpl]';
        detectedPplDEBUGGING = indexVecForAllPpl(canBeDetectedDEBUGGING & detectedDEBUGGING);
        diagnosedPool(timePeriod,1) = size(activeAllPpl,1);  %not on trt (active people only if ssspec = 1)
        diagnosedPool(timePeriod,2) = size(coveredAllPpl_maySeek,1);  %also covered by RNTCP
        diagnosedPool(timePeriod,3) = size(coveredAllPpl,1);  %willing to seek RNTCP help
        diagnosedPool(timePeriod,4) = sum(canBeDetectedDEBUGGING);  %number of people showing up for test calibrated
        diagnosedPool(timePeriod,5) = sum(detectedDEBUGGING);  %would have tested pos
        diagnosedPool(timePeriod,6) = size(detectedPplDEBUGGING,1); %showed up for test and tested pos
                
        detected = rand(numCurrentPpl,1) < detectionProb;
        indexVecForAllPpl = [1:numCurrentPpl]';
        detectedPpl = indexVecForAllPpl(detected); %subset active people who are detected by test
        detectedPpl_TB = detectedPpl;
        detectedPpl_noTBageSex = [];
        totalScreened = (size(detectedPpl,1)*(1/TBparams.SSsensit));
        
        notification_TBmac(timePeriod,1) = sum(notifications(detectedPpl,1) == 0);
        detectedPpl_over15 = subsetterFuncInequality(stateMat,detectedPpl,Page,16,109);         
        notification_TBmac(timePeriod,5) = sum(notifications(detectedPpl_over15,1) == 0);
        notifications(detectedPpl,1) = 1;  %these people have been detected

        if TBparams.SSspec < 1
            %get the noTB people who are detected because spec < 1
            detectedPpl_TB = subsetterFunc(stateMat,detectedPpl,Phealth,[3,4]); %active people detected by test
            detectedPpl_noTB = subsetterFunc(stateMat,detectedPpl,Phealth,[0,1,2]); %inactive ppl detected by test
            numTurnedAwwwayPerDiagnosedVar = TBparams.numTurnedAwayPerDiagnosed;
            for ageGrp = 1:1:102
                detectedPpl_inAgeNoTB = subsetterFunc(stateMat,detectedPpl_noTB,Page,[ageGrp]); %noTB but diagnosed in age group
                detectedPpl_inAgeTB = subsetterFunc(stateMat,detectedPpl_TB,Page,[ageGrp]); %TB, diagnosed in age group
                for sexGrp = 1:1:2
                    detectedPpl_inSexNoTB = subsetterFunc(stateMat,detectedPpl_inAgeNoTB,Psex,[sexGrp]); %and same sex
                    detectedPpl_inSexTB = subsetterFunc(stateMat,detectedPpl_inAgeTB,Psex,[sexGrp]); %and same sex
                    probNoTBDiag = (numTurnedAwwwayPerDiagnosedVar*size(detectedPpl_inSexTB,1))/size(detectedPpl_inSexNoTB,1);
                    ADetectedPpl_noTB = rand(size(detectedPpl_inSexNoTB,1),1) <= probNoTBDiag;
                    detectedPpl_noTBageSex = [detectedPpl_noTBageSex; detectedPpl_inSexNoTB(ADetectedPpl_noTB)];
                end
            end
            detectedPpl_noTBageSex = detectedPpl_noTBageSex(detectedPpl_noTBageSex~=0);  %detectedPpl_noTBageSex now has no zeros
            totalScreened = (size(detectedPpl_TB,1)*(1/TBparams.SSsensit)) + (size(detectedPpl_noTBageSex,1)*(1/(1-TBparams.SSspec)));
            detectedPpl = unique(sort([detectedPpl_TB;detectedPpl_noTBageSex]));
        end
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + ( totalScreened*TBparams.SStestCost );  %every screened person has SS cost.  TBparams.numTurnedAwayPerDiagnosed = For every diagnosed person, there are this many people who seek treatment but turned away because test neg
        
        if (TBMAC4_activeCaseFinding == 1 || TBMAC5_preventativeTherapyForLat == 1)  %active case finding scenario
            startACFtime = 1788;  stopACFtime = 1849; %1788 is 2015, 1800 is 2016, so 2018 is timePeriod 1825, and 2020 is timePeriod1849, 2025 is timePeriod1908
            if Tbmac4_adv == 1
               stopACFtime = 1909;
            end
            if (timePeriod >= startACFtime && timePeriod < stopACFtime && mod((timePeriod-1788),TBmac4_freq)==0)  %ends in 2025
                %actively screen some people and put them on treatment
                aliveNow = subsetterFunc(stateMat,[-1],Phealth,[0,1,2,3,4]); % not on treatment
                noTrt = subsetterFunc(stateMat,aliveNow,Ptreatment,0); % not on treatment
                screened = ( rand(size(noTrt,1),1)) < ACFpercOfUntrted;  %screen this proportion of the untreated population every year
                numScreened(timePeriod,1) = sum(screened);
                screenedHasTB = subsetterFunc(stateMat,noTrt(screened),Phealth,[3,4]); % has TB
                if TBMAC5_preventativeTherapyForLat == 1
                    screenedNoTb = subsetterFunc(stateMat,noTrt(screened),Phealth,[0]);    %is healthy
                else
                    screenedNoTb = setdiff(noTrt(screened), screenedHasTB);    %does not have TB
                end
                ACFfound_TB = (rand(size(screenedHasTB,1),1)) < TBmac4_sens; %TBparams.geneXsensDS;  %detected by geneX
                ACFfound_noTB = ( rand(size(screenedNoTb,1),1)) < TBmac4_spec; %TBparams.geneXspecDS;  %erroneously detected by geneX.  make spec perfect
                detectedPpl_TB = sort([detectedPpl_TB; screenedHasTB(ACFfound_TB)]);
                detectedPpl_noTBageSex = sort([detectedPpl_noTBageSex; screenedNoTb(ACFfound_noTB)]);
                numScreened(timePeriod,2) = sum(ACFfound_TB);
                
                if TBMAC5_preventativeTherapyForLat == 1
                    numScreened(timePeriod,3) = numScreened(timePeriod,1) - size(screenedHasTB,1);  %detected and start LTBI trt
                    latentScreened = subsetterFunc(stateMat,noTrt(screened),Phealth,[1,2]); % screened pool and has latent TB
                    fracScreened = 1;
                    ACFfound_LTBI = (rand(size(latentScreened,1),1)) < (fracScreened*TBmac4_sens*0.82*0.8*TBparams.TBmacMultiplier);  %detected by latent TB test (TST) x proportion completed (82%) * % protection offered (80%). 
                    curable = subsetterFunc(stateMat,latentScreened(ACFfound_LTBI),Phealth,[1]); % screened pool and has latent TB
                    stateMat(curable, Phealth) = 0; %cured
                    numScreened(timePeriod,4) = size(latentScreened,1)*fracScreened*TBmac4_sens;  %detected and start LTBI trt
                end
            end
        end
        
        %offset detected people: wait a month if just activated
        detectedPpl_notTooEarly = subsetterFuncInequality(stateMat,detectedPpl_TB,PactTime,0,timePeriod-1); %activated before this time period, eligible for trtmt
        detectedPpl_tooEarly = setdiff(detectedPpl_TB,detectedPpl_notTooEarly);
        waitingForTrt = subsetterFunc(stateMat,[-1],PtrtmtQ,1); %are in queue for treatment
        stateMat(waitingForTrt, PtrtmtQ) = 0; %reset their flag
        detectedPpl = sort(unique([detectedPpl_notTooEarly;waitingForTrt;detectedPpl_noTBageSex;privateToRNTCP]));
        privateToRNTCP = [];
        stateMat(detectedPpl_tooEarly, PtrtmtQ) = 1;  %too-early people wait a month if detected
        stateMat(detectedPpl, PtrtmtQ_private) = 0;
        
        %for the privateClinicInsteadOfRNTCP scenario
        if TBparams.privateClinicInsteadOfRNTCP == 1
            detectedInRNTCPPpl = detectedPpl;
            detectedPpl = [];
        end
        
        if TBparams.DSTinsteadOfSS  == 1
            %want to do DST right away
            stateMat(detectedPpl,PmdrTesting)       = timePeriod;
            stateMat(detectedPpl,PmdrTestingStatus) = stateMat(detectedPpl, Phealth);
            %             mdrTestPpl(timePeriod, 13) = size(subsetterFunc(stateMat, detectedPpl, Phealth, [0,1,2]),1);
            %             mdrTestPpl(timePeriod, 14) = size(subsetterFunc(stateMat, detectedPpl, Phealth, [3]),1);
            %             mdrTestPpl(timePeriod, 15) = size(subsetterFunc(stateMat, detectedPpl, Phealth, [4]),1);
        end
        
        if TBMAC3_XpertReplacesSmear == 1
            if timePeriod > 1800
                coveredDetected_TBmac = rand(size(detectedPpl,1),1) < GeneXcoverageMat(min(1848-1800, timePeriod-1800));  %undergo MDR diagnosis right away
                stateMat(detectedPpl(coveredDetected_TBmac),PmdrTesting)       = timePeriod;
                stateMat(detectedPpl(coveredDetected_TBmac),PmdrTestingStatus) = stateMat(detectedPpl(coveredDetected_TBmac), Phealth);
            end
        end
        
        %count the diagnosed ppl to make some statistics
        diagnosisCounterMat(timePeriod, 1) = sum(updating);          % total people currently being simulated
        diagnosisCounterMat(timePeriod, 2) = size(activeTBppl_all,1);% total
        diagnosisCounterMat(timePeriod, 3) = size(activeAllPpl, 1);   % actTB ppl not on treatment
        diagnosisCounterMat(timePeriod, 4) = size(detectedPpl,1);          % number ppl put on treatment
        activeTB_putOnTrt = subsetterFunc(stateMat,detectedPpl,Phealth,[3,4]); %active people
        timeSinceAct = timePeriod - stateMat(activeTB_putOnTrt, PactTime);
        diagnosisCounterMat(timePeriod, 8) = mean(timeSinceAct);     % average num months since activation
        ignoredTBtime = [ignoredTBtime;timeSinceAct];
        %DEBUGGING
        %diagnosedAges(timePeriod,1:21) = hist( stateMat(detectedPpl, Page),[0:5:100]);
        
        %Retreatment go to Cat II
        retreatment = stateMat(detectedPpl,PpastTreatment ) > 0;
        stateMat(detectedPpl( retreatment), Ptreatment) = 2; % Retreatment go to CATII if they are detected
        diagnosisCounterMat(timePeriod, 7) = size(detectedPpl( retreatment),1);
        diagnosisCounterMat(timePeriod, 10) = sum(stateMat(detectedPpl(retreatment),Phealth) == 4); %mdr ppl put in cat II
        diagnosisCounterMat(timePeriod, 11) = sum(stateMat(detectedPpl(retreatment),PpastTreatment) == 3); %relapse ppl put in cat II
        diagnosisCounterMat(timePeriod, 12) = sum(stateMat(detectedPpl(retreatment),PpastTreatment) == 3 & stateMat(detectedPpl(retreatment),Phealth) == 4); %mdr relapse ppl put in cat II
        
        %get DST for all cat II  (TBMAC modification to base case)
        if timePeriod >= 1800
            fracRetrt = rand(size(detectedPpl,1),1) <= min(1, ( double(timePeriod-1800) )/60);
            stateMat(detectedPpl(fracRetrt),PmdrTesting)       = timePeriod;
            stateMat(detectedPpl(fracRetrt),PmdrTestingStatus) = stateMat(detectedPpl(fracRetrt), Phealth);
        end
        
        %nonretreatment
        stateMat(detectedPpl(~retreatment), Ptreatment) = 1; % Non-retreaters go to CatI if they are detected
        diagnosisCounterMat(timePeriod, 6) = size(detectedPpl(~retreatment),1);
        diagnosisCounterMat(timePeriod, 9) = sum(stateMat(detectedPpl(~retreatment),Phealth) == 4); %mdr ppl put in cat I
        
        %write time to diagnosis for first time treatment people
        timeSinceAct = timePeriod - stateMat(detectedPpl(~retreatment), PactTime);
        diagnosisTimeMat(timePeriod, 1) = mean(timeSinceAct); % average num months since activation overall for CATI
        if isempty(timeSinceAct) == 0
            diagnosisTimeMat(timePeriod, 2) = median(single(timeSinceAct));
        end
        %         COMMENTED OUT FOR SPEED
        %         diagnosisTimeDist(timePeriod,1:20) = histc(timeSinceAct,[1:1:20])';
        %         counter = 3;
        %         desiredAge = [0:10:90];
        %         for j = 1:size(desiredAge,2)
        %             detectedUnder40 = subsetterFuncInequality(stateMat,detectedPpl(~retreatment),Page,desiredAge(j),desiredAge(j)+4); %under40
        %             for i = 1:2
        %                 desiredSex = [1 , 2];
        %                 detectedUnder40Male = subsetterFunc(stateMat,detectedUnder40,Psex,desiredSex(i)); %males
        %                 timeSinceAct = timePeriod - stateMat(detectedUnder40Male, PactTime);
        %                 diagnosisTimeMat(timePeriod, counter) = mean(timeSinceAct); % average num months since activation
        %                 if isempty(timeSinceAct) == 0
        %                     diagnosisTimeMat(timePeriod, counter+1) = median(single(timeSinceAct)); %median time
        %                 end
        %                 diagnosisTimeMat(timePeriod, counter+2) = size(detectedUnder40Male,1); %num in grp
        %                 counter = counter + 3;
        %             end
        %         end
        
        %treatment flag for everyone detected (all some treatment)
        stateMat(detectedPpl, PpastTreatment) = 1; %treatment flag;
        
        %DOTS PLUS DIAGNOSIS
        if (timePeriod > DotsPlusStartTime) %one year after rntcp full coverage
            dotsPlusCoverage = TBparams.DotsPlusFullCoverage;
            %change to less coverage if not fully implemented (prior to 2012)
            if (timePeriod < DotsPlusStartTime + TBparams.DotsPlusRampUpPeriod )
                dotsPlusCoverage = TBparams.DotsPlusRampUpSequence(timePeriod-DotsPlusStartTime,1)*0.13;  %modified for TBmac
            end
            dotsPlusCoverage
            MDRtestingDone    = subsetterFuncInequality(stateMat,[-1],PmdrTesting ,1,timePeriod-TBparams.DSTperiodLength);  %tested for MDR over 6 months ago
            %MDRtestingDone     = subsetterFunc(stateMat,MDRtestingDone1,Ptreatment,[0,1,2]);  %not already on treatment
            logical_MDRcovered = rand(size(MDRtestingDone,1),1) <= dotsPlusCoverage; %need to check for coverage
            MDRcoveredAll      = MDRtestingDone(logical_MDRcovered);
            MDRcovered         = subsetterFunc(stateMat, MDRcoveredAll, Phealth, [4]); %has MDR: LJ test gold standard, assume gets all MDR statuses correct
            
            notification_TBmac(timePeriod,2) =  size(MDRtestingDone,1);

            %if DST test not 100% specific, need to do:
            MDRtestingDoneNoMDR = subsetterFunc(stateMat, MDRcoveredAll, Phealth, [0,1,2,3]);
            testedPosForMDR_noMDR = rand(size(MDRtestingDoneNoMDR,1),1) <= (1-TBparams.DSTspec);  %got all these people wrong.  CHECK IF 1- or not.
            
            %if DST test not 100% sensitive, need to do:
            testedPosForMDR_mdr = rand(size(MDRcovered,1),1) <= TBparams.DSTsensit;
            
            additionalDiagnosed_tbmac = [];
            if TBMAC3_XpertReplacesSmear == 1
                if timePeriod > 1800
                    coveredGeneX_TBmac = rand(size(MDRcovered,1),1) < GeneXcoverageMat(min(1848-1800, timePeriod-1800));  %get geneX test
                    geneXcovered = MDRcovered(coveredGeneX_TBmac);
                    testedPosMDR_tbMac = rand(size(geneXcovered,1),1) <= TBparams.geneXsensMDR;
                    additionalDiagnosed_tbmac = geneXcovered(testedPosMDR_tbMac);
                end
            end
            
            %put people on treatment
            enterCatIV = unique(sort([MDRtestingDoneNoMDR(testedPosForMDR_noMDR);MDRcovered(testedPosForMDR_mdr);additionalDiagnosed_tbmac]));
            stateMat(enterCatIV,Ptreatment) = 4;
            stateMat(enterCatIV,PtrtmtCounter) = 0;
            allTested = unique(sort([MDRtestingDoneNoMDR;MDRcovered]));
            stateMat(allTested, PmdrTesting) = 0; %reset to default
            stateMat(allTested, PmdrTestingStatus) = 0;
            
            notification_TBmac(timePeriod,3) =  sum(notifications(enterCatIV,1) == 0);
            notifications(enterCatIV,1) = 1;  %these people have been detected
            
            %             mdrTestPpl(timePeriod,1) = size(MDRtestingDoneNoMDR,1);
            %             mdrTestPpl(timePeriod,2) = sum(testedPosForMDR_noMDR);
            %             mdrTestPpl(timePeriod,3) = size(MDRcovered,1);
            %             mdrTestPpl(timePeriod,4) = sum(testedPosForMDR_mdr);
            %
            %DEBUGGING APR 2014
            %             MDRspecPool_healthy = subsetterFunc(stateMat, MDRtestingDoneNoMDR, Phealth, [0,1,2]);
            %             MDRspecPool_healthyNoTrt = subsetterFunc(stateMat, MDRspecPool_healthy, Ptreatment, [0]);
            %             MDRspecPool_healthyDOTS = subsetterFunc(stateMat, MDRspecPool_healthy, Ptreatment, [1,2]);
            %             MDRspecPool_DS = subsetterFunc(stateMat, MDRtestingDoneNoMDR, Phealth, [3]);
            %             MDRspecPool_DSNoTrt = subsetterFunc(stateMat, MDRspecPool_DS, Ptreatment, [0]);
            %             MDRspecPool_DSDOTS = subsetterFunc(stateMat, MDRspecPool_DS, Ptreatment, [1,2]);
            %
            %             mdrTestPpl(timePeriod,5) = size(MDRspecPool_healthy,1);
            %             mdrTestPpl(timePeriod,6) = size(MDRspecPool_healthyNoTrt,1);
            %             mdrTestPpl(timePeriod,7) = size(MDRspecPool_healthyDOTS,1);
            %             mdrTestPpl(timePeriod,8) = mean(double(stateMat(MDRspecPool_healthyDOTS,PtrtmtCounter)));
            %             mdrTestPpl(timePeriod,9) = size(MDRspecPool_DS,1);
            %             mdrTestPpl(timePeriod,10) = size(MDRspecPool_DSNoTrt,1);
            %             mdrTestPpl(timePeriod,11) = size(MDRspecPool_DSDOTS,1);
            %             mdrTestPpl(timePeriod,12) = mean(double(stateMat(MDRspecPool_DSDOTS,PtrtmtCounter)));
            %
            
            
            %count the diagnosed ppl to make some statistics
            diagnosisCounterMat(timePeriod, 5) = size(enterCatIV,1); %num MDR ppl put on treatment
        end
    end  %end "if (timePeriod >= burnInTime)"
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'normal diagnosis';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%%%%% PRIVATE CLINIC DIAGNOSIS AND MDR %%%%%%%%%%%%%%%%%%
    MdrCasesEvolveAct_private = 0;
    detectedPpl_private = [];
    if (timePeriod >= burnInTime) %after burn in
        %get eligible people
        if TBparams.SSspecPrivate < 1
            activeTBppl_notRNTCPb4 = subsetterFunc(stateMat,[-1],PpastTreatment,[0]); %never had RNTCP treatment %not just active TB ppl if spec less than one
        else
            activeTBppl_all = subsetterFunc(stateMat,-1,Phealth,[3,4]); %active people
            activeTBppl_notRNTCPb4 = subsetterFunc(stateMat,activeTBppl_all,PpastTreatment,[0]); %never had RNTCP treatment
        end
        activeTBppl_elig = subsetterFuncInequality(stateMat,activeTBppl_notRNTCPb4,PprivTrtCount,0,TBparams.maxPrivateTrtCount); %active people
        activeTBppl     = subsetterFunc(stateMat,activeTBppl_elig,Ptreatment,[0]); %active people not on treatment
        
        %identify TB knowledge for each person
        %same as above.  knowledgeAgeBrac = ageBracketMaker(TBparams.uptakeAgeBracs, stateMat(1:numCurrentPpl,Page));
        % just use urbanKnowl from above.  treatmentKnowl = [TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,1), TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,2)];
        knowledgeAgeBrac = ageBracketMaker(TBparams.uptakeAgeBracs, stateMat(1:numCurrentPpl,Page));
        urbanKnowl = [TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,1), TBparams.uptakeUrbanKnowledge(knowledgeAgeBrac,2)];
        TBknowledge = -1*ones(numCurrentPpl,1);
        for ii = 1:size(activeTBppl,1)
            TBknowledge(activeTBppl(ii,1),1) = urbanKnowl(activeTBppl(ii,1),stateMat(activeTBppl(ii,1),Psex) );
        end
        
        %get the proper detection prob
        sensitUnique = ones(numCurrentPpl, 1)*(1-TBparams.SSspecPrivate);
        activeTBppl_all = subsetterFunc(stateMat,[-1],Phealth,[3,4]); %active people
        sensitUnique(activeTBppl_all) = TBparams.SSsensitPriv;
        
        %display out
        if timePeriod == burnInTime + 120
            disp('this is the private clinic prob')
            TBparams.prTestedGivenTB*TBparams.seekPrivate*TBparams.SSsensitPriv.*TBparams.uptakeUrbanKnowledge
            TBparams.prTestedGivenTB*TBparams.seekPrivate*(1-TBparams.SSspecPrivate).*TBparams.uptakeUrbanKnowledge
        end
        
        %diagnose people
        detectionProb = TBparams.prTestedGivenTB*TBparams.seekPrivate*sensitUnique.*TBknowledge;
        if any(detectionProb > 1)
            detectionProb(detectionProb > 1) = 1;
        end
        if any(detectionProb < 0)
            detectionProb(detectionProb < 0) = 0;
        end
        detected = rand(numCurrentPpl,1) < detectionProb;
        indexVecForAllPpl = [1:numCurrentPpl]';
        detectedPpl = indexVecForAllPpl(detected); %subset active people who are detected by test
        detectedPpl_TB = detectedPpl;
        totalScreened = (size(detectedPpl,1)*(1/TBparams.SSsensitPriv));
        
        notification_TBmac(timePeriod,4) =   sum(notifications(detectedPpl,1) == 0);
        detectedPpl_over15 = subsetterFuncInequality(stateMat,detectedPpl,Page,16,109);
        notification_TBmac(timePeriod,6) = sum(notifications(detectedPpl_over15,1) == 0);
        notifications(detectedPpl,1) = 1;  %these people have been detected
        
        notification_TBmac(timePeriod,7) = sum(privateTrtNaive(detectedPpl,1) == 0);  %count how many trt naive diagnosed (not notified)
        privateTrtNaive(detectedPpl,1) = 1;  %mark these ppl as been in private trt
        
        if TBparams.SSspecPrivate < 1
            %get the noTB people who are detected because spec < 1
            detectedPpl_TB = subsetterFunc(stateMat,detectedPpl,Phealth,[3,4]); %active people detected by test
            detectedPpl_noTB = subsetterFunc(stateMat,detectedPpl,Phealth,[0,1,2]); %inactive ppl detected by test
            detectedPpl_noTBageSex = [];
            numTurnedAwwwayPerDiagnosedVar = TBparams.numTurnedAwayPerDiagnosed;
            parfor ageGrp = 1:1:102
                detectedPpl_inAgeNoTB = subsetterFunc(stateMat,detectedPpl_noTB,Page,[ageGrp]); %noTB but diagnosed in age group
                detectedPpl_inAgeTB = subsetterFunc(stateMat,detectedPpl_TB,Page,[ageGrp]); %TB, diagnosed in age group
                for sexGrp = 1:1:2
                    detectedPpl_inSexNoTB = subsetterFunc(stateMat,detectedPpl_inAgeNoTB,Psex,[sexGrp]); %and same sex
                    detectedPpl_inSexTB = subsetterFunc(stateMat,detectedPpl_inAgeTB,Psex,[sexGrp]); %and same sex
                    probNoTBDiag = (numTurnedAwwwayPerDiagnosedVar*size(detectedPpl_inSexTB,1))/size(detectedPpl_inSexNoTB,1);
                    ADetectedPpl_noTB = rand(size(detectedPpl_inSexNoTB,1),1) <= probNoTBDiag;
                    detectedPpl_noTBageSex = [detectedPpl_noTBageSex; detectedPpl_inSexNoTB(ADetectedPpl_noTB)];
                end
            end
            detectedPpl_noTBageSex = detectedPpl_noTBageSex(detectedPpl_noTBageSex~=0);  %detectedPpl_noTBageSex now has no zeros
            totalScreened = (size(detectedPpl_TB,1)*(1/TBparams.SSsensit)) + (size(detectedPpl_noTBageSex,1)*(1/(1-TBparams.SSspec)));
            detectedPpl = unique(sort([detectedPpl_TB;detectedPpl_noTBageSex]));
        end
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + ( totalScreened*TBparams.PrivateTestCost);  %every screened person has SS cost.  TBparams.numTurnedAwayPerDiagnosed = For every diagnosed person, there are this many people who seek treatment but turned away because test neg
        
        %offset detected people: wait a month if just activated
        detectedPpl_notTooEarly = subsetterFuncInequality(stateMat,detectedPpl_TB,PactTime,0,timePeriod-1); %activated before this time period, eligible for trtmt
        detectedPpl_tooEarly = setdiff(detectedPpl_TB,detectedPpl_notTooEarly);
        waitingForTrt = subsetterFunc(stateMat,[-1],PtrtmtQ_private,1); %are in queue for treatment
        stateMat(waitingForTrt, PtrtmtQ_private) = 0; %reset their flag
        detectedPpl = sort(unique([detectedPpl_notTooEarly;waitingForTrt;detectedPpl_noTBageSex]));
        stateMat(detectedPpl_tooEarly, PtrtmtQ_private) = 1;  %too-early people wait a month if detected
        stateMat(detectedPpl, PprivTrtCount) = stateMat(detectedPpl, PprivTrtCount) + 1;  %onPrivate trt
        detectedPpl_private = detectedPpl;
        
        %count ppl being put on trt
        timeSinceHadTB = timePeriod - stateMat(detectedPpl, PactTime);  %count these people
        ignoredTBtime = [ignoredTBtime;timeSinceHadTB];  %write down duration of ignored TB
        
        %for the privateClinicInsteadOfRNTCP scenario
        if TBparams.privateClinicInsteadOfRNTCP == 1
            detectedPpl = sort(unique([detectedPpl_notTooEarly;waitingForTrt;detectedInRNTCPPpl]));
            detectedInRNTCPPpl = [];
        end
        
        %for the private-public mix scenario, where private refers to public
        if TBparams.privateReferToRNTCP ~= 0
            referredToRNTCP = rand(length(detectedPpl),1) <= (TBparams.privateReferToRNTCP * TBparams.PPMeffectiveness );
            costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) +((length(detectedPpl)*TBparams.privateReferToRNTCP)*TBparams.perPatientPPMcost);  %PPM cost if try to refer patients
            privateToRNTCP = detectedPpl(referredToRNTCP);
            detectedPpl = detectedPpl(~referredToRNTCP);
            numRefPpl(timePeriod-burnInTime+1,1) =  size(detectedPpl_private,1);  %total entering private this month
            numRefPpl(timePeriod-burnInTime+1,2) = size(detectedPpl,1);  %still in private treatment
            numRefPpl(timePeriod-burnInTime+1,3) = size(privateToRNTCP,1);
        end
        
        if TBMAC1b_increaseTrt == 1  %ramp up 80% ppm over 5 years (120 months) starting in 2015
            if timePeriod > 1800
                referredToRNTCP = rand(length(detectedPpl),1) <= TBparams.PPMeffectivenessVal(min(1884-1800, timePeriod-1800)) ;
                privateToRNTCP = detectedPpl(referredToRNTCP);
                detectedPpl = detectedPpl(~referredToRNTCP);
                numRefPpl(timePeriod-burnInTime+1,1) =  size(detectedPpl_private,1);  %total entering private this month
                numRefPpl(timePeriod-burnInTime+1,2) = size(detectedPpl,1);  %still in private treatment
                numRefPpl(timePeriod-burnInTime+1,3) = size(privateToRNTCP,1);
            end
        end
        
        %count the diagnosed ppl to make some statistics
        privateCounterMat(timePeriod, 1) = sum(updating);          % total people currently being simulated
        privateCounterMat(timePeriod, 2) = size(activeTBppl_all,1);% total active TB ppl
        privateCounterMat(timePeriod, 3) = size(activeTBppl, 1);   % actTB ppl not on treatment
        privateCounterMat(timePeriod, 4) = size(detectedPpl, 1);          % number ppl put on treatment
        timeSinceAct = timePeriod - stateMat(detectedPpl, PactTime);
        privateCounterMat(timePeriod, 5) = mean(timeSinceAct);     % average num months since activation
        
        %for scenarios where can get cured by private
        if TBparams.privateClinicCure ~= 0
            DStbPpl = subsetterFunc(stateMat,detectedPpl,Phealth,[3]); %active DS people
            privateCured = rand(size(DStbPpl,1),1) < TBparams.privateClinicCure;
            stateMat(DStbPpl(privateCured),Phealth)        = 0;
            detectedPpl = detectedPpl(~privateCured);
        end
        
        %in this month on private treatment you may get mdr
        costsMat(timePeriod-burnInTime+1,1) = costsMat(timePeriod-burnInTime+1,1) + (size(detectedPpl,1)*(TBparams.monthlyPrivClinicCost-TBparams.SStestCost));  %private treatment ppl incur priv treatment cost
        gotMDR = rand(size(detectedPpl,1),1) < TBparams.privateMDRprob;
        stateMat(detectedPpl(gotMDR), Phealth) = 4;  %get MDR
        MdrCasesEvolveAct_private = sum(gotMDR);
        numPrivateClinics(timePeriod,1:8) = histc(stateMat(1:numCurrentPpl, PprivTrtCount),[1:1:8])';
    end  %end after burn in time
    
    %%%%%%%%%% SELF CURE
    if TBparams.probSelfCure > 0  %monthly probability
        %monthly self cure
        TBppl_notPriv = subsetterFunc(stateMat,setdiff([1:numCurrentPpl]',detectedPpl_private),Phealth,[3,4]);  %those not on private trt and TB
        %         allActPplForCure = subsetterFunc(stateMat,TBppl_notPriv,Ptreatment,[0,1]);  %untreated or RNTCP catI
        allActPplForCure = subsetterFunc(stateMat,TBppl_notPriv,Ptreatment,[0]);  %untreated
        
        %some fraction of them go back to healthy
        selfCured = rand(length(allActPplForCure),1) <= TBparams.probSelfCure;
        stateMat(allActPplForCure(selfCured),Phealth) = 0; % Set their health state to healthy
        %       trtCountMat(timePeriod,2) = sum(selfCured);  %record the num selfCured
        
        timeSinceHadTB = timePeriod - stateMat(allActPplForCure(selfCured), PactTime);  %count these people
        ignoredTBtime = [ignoredTBtime;timeSinceHadTB];  %write down duration of ignored TB
    end
    
    
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'private treatment';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    %%%%%%%%RECORD THE MOST RECENT STATE%%%%%%%%%%%
    
    if (timePeriod >= burnInTime)
        ignoredTBmat(timePeriod, 1) = mean(ignoredTBtime);  %record the mean duration of ignored TB, in months
        
        %get the number of people for QALYS
        healthyLat = subsetterFunc(stateMat,[-1],Phealth,[0,1,2]);
        healthyLatDOTS = size(subsetterFunc(stateMat,healthyLat,Ptreatment,[1,2]),1);  %get QALY = 0.851699
        healthyLatMDRTrt = size(subsetterFunc(stateMat,healthyLat,Ptreatment,[4]),1); %get QALY = 0.806699
        healthyLatNoTrt = subsetterFunc(stateMat,healthyLat,Ptreatment,[0]);
        healthyLatNeverTrt = size(subsetterFunc(stateMat,healthyLatNoTrt,PpastTreatment,[0]),1);  %get QALY = 1
        healthyLatPastTrt = size(healthyLatNoTrt,1) - healthyLatNeverTrt; %get QALY = 0.941699
        
        sickPpl = subsetterFunc(stateMat,[-1],Phealth,[3,4]);
        
        DSppl = subsetterFunc(stateMat,sickPpl,Phealth,[3]);
        DSpplNoTrt = size(subsetterFunc(stateMat,DSppl,Ptreatment,[0]),1);  %get 0.663
        DSpplInMDR = size(subsetterFunc(stateMat,DSppl,Ptreatment,[4]),1);  %get 0.753
        DSpplInDOTS = size(DSppl,1) - DSpplNoTrt - DSpplInMDR;  %in cat I and II.  get 0.843
        
        MDRppl = subsetterFunc(stateMat,sickPpl,Phealth,[4]);
        MDRpplInMDR = size(subsetterFunc(stateMat,MDRppl,Ptreatment,[4]),1);  %get 0.753
        MDRpplNotInMDR = size(MDRppl,1) - MDRpplInMDR ; %get 0.663
        
        numQalyPpl(1:9,timePeriod-burnInTime+1) = [healthyLatNeverTrt;healthyLatPastTrt;healthyLatDOTS;healthyLatMDRTrt;DSpplNoTrt;DSpplInDOTS;DSpplInMDR;MDRpplNotInMDR;MDRpplInMDR];
        numQalyPpl(19,timePeriod-burnInTime+1) = numCurrentPpl;  % sum(updating);  %numCurrentPpl
    end
    
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'QALYs';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if LEbuilder == 1
        LEnumQalyPpl(:,timePeriod-origTotPeriods) = numQalyPpl(:,timePeriod-burnInTime+1);
        LEcostsMat(timePeriod-origTotPeriods,1) = costsMat(timePeriod-burnInTime+1,1);
    elseif durationYrs ~= 230
        if ( timePeriod > startRecordingTime && mod(timePeriod-startRecordingTime,recordUnit) == 0)  %record only if the Nth year
            %recNo has all the info.  Order, from innermost: PurbanRur,PpastTreatment, smoking,
            %age, sex, health, treatment
            doubleStateMat = double(stateMat);
            recNo(:,((timePeriod-startRecordingTime)/recordUnit)) =  (((((((( ( ( doubleStateMat(:,PurbanRur))*3) + ...
                doubleStateMat(:,Psmoking))*101) + ...
                doubleStateMat(:,Page))*3) + doubleStateMat(:,Psex))*8) + ...
                doubleStateMat(:,Phealth))*6) + doubleStateMat(:,Ptreatment);
        end
        
        %treatment and health states recording COMMENTED OUT FOR SPEED
        %     if mod(timePeriod, 10)==0
        %         %make health state count for people in each treatment cat
        %         peoples = subsetterFunc(stateMat,[-1],Ptreatment,[0]);
        %         NoTreatHealthStates(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,[-1],Ptreatment,[1]);
        %         CatIHealthStates(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,[-1],Ptreatment,[2]);
        %         CatIIHealthStates(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,[-1],Ptreatment,[4]);
        %         CatIVHealthStates(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %
        %         %updating only
        %         peoples = subsetterFunc(stateMat,find(updating==1),Ptreatment,[0]);
        %         NoTreatHealthStates_updating(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,find(updating==1),Ptreatment,[1]);
        %         CatIHealthStates_updating(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,find(updating==1),Ptreatment,[2]);
        %         CatIIHealthStates_updating(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %         peoples = subsetterFunc(stateMat,find(updating==1),Ptreatment,[4]);
        %         CatIVHealthStates_updating(timePeriod/10,:) = histc(stateMat(peoples, Phealth), [0:1:7]);
        %     end
        %end treatment and health states recording
        
        %incidence graph numbers
        if timePeriod < burnInTime
            MdrCasesEvolvedAct = 0;
        end
        
        %incidenceMat(timePeriod,1) = activaCounterMat(timePeriod, 3)+activaCounterMat(timePeriod, 4)+activaCounterMat(timePeriod, 5)+activaCounterMat(timePeriod, 6) + MdrCasesEvolvedAct + MdrCasesEvolveAct_private;
        %not counting DS to MDR as a separate incidence case
        incidenceMat(timePeriod,1) = activaCounterMat(timePeriod, 3)+activaCounterMat(timePeriod, 4)+activaCounterMat(timePeriod, 5)+activaCounterMat(timePeriod, 6);
        %    incidenceMat(timePeriod,2) = incidenceMat(timePeriod,2) + sum(MDRtransmitted);  %not doing right now
        incidenceMat(timePeriod,3) = sum(  stateMat(1:numCurrentPpl, Phealth) ~= 5 & stateMat(1:numCurrentPpl, Phealth) ~= 6   ); %also get the number of alive for graph
        
        %write down the disease cases and deaths for averting stuff
        if (timePeriod >= burnInTime)
            %monthly actTB prev by category
            monthlyActOutcomes(timePeriod-burnInTime+1, 1) = sum(stateMat(:,Phealth) == 0 & stateMat(:,Ptreatment) == 0 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 2) = sum(stateMat(:,Phealth) == 1 & stateMat(:,Ptreatment) == 0 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 3) = sum(stateMat(:,Phealth) == 2 & stateMat(:,Ptreatment) == 0 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 4) = sum(stateMat(:,Phealth) == 3 & stateMat(:,Ptreatment) == 0 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 5) = sum(stateMat(:,Phealth) == 3 & stateMat(:,Ptreatment) == 1 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 6) = sum(stateMat(:,Phealth) == 3 & stateMat(:,Ptreatment) == 2 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 7) = sum(stateMat(:,Phealth) == 3 & stateMat(:,Ptreatment) == 4 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 8) = sum(stateMat(:,Phealth) == 4 & stateMat(:,Ptreatment) == 0 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 9) = sum(stateMat(:,Phealth) == 4 & stateMat(:,Ptreatment) == 1 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 10) = sum(stateMat(:,Phealth) == 4 & stateMat(:,Ptreatment) == 2 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 11) = sum(stateMat(:,Phealth) == 4 & stateMat(:,Ptreatment) == 4 );
            monthlyActOutcomes(timePeriod-burnInTime+1, 12) =  sum(updating);
            
            %now only counting active MDR incidence
            %MdrMethod(timePeriod-burnInTime+1,1) = MdrCasesEvolvedAct + MdrCasesEvolvedLat;
            %MdrMethod(timePeriod-burnInTime+1,2) = MdrCasesTransmitted;
            MdrMethod(timePeriod-burnInTime+1,1) = MdrCasesEvolvedAct + MdrCasesEvolvedLatNowAct;
            MdrMethod(timePeriod-burnInTime+1,2) = MdrCasesActivated - MdrCasesEvolvedLatNowAct;
            MdrMethod(timePeriod-burnInTime+1,3) =  sum(updating);
            MdrMethod(timePeriod-burnInTime+1,4) = MdrCasesEvolveAct_private;
            MdrActivation(timePeriod-burnInTime+1,1) = MdrCasesActivated;
            
            diseaseCasesMat(timePeriod-burnInTime+1, 1) = DSCasesTransmitted; %latent DS infection incidence
            diseaseCasesMat(timePeriod-burnInTime+1, 2) = MdrCasesEvolvedLat + MdrCasesTransmitted; %latent MDR infection incidence
            diseaseCasesMat(timePeriod-burnInTime+1, 3) = sum(activaCounterMat(timePeriod, 3:6)) - MdrCasesActivated; %active DS disease incidence;
            diseaseCasesMat(timePeriod-burnInTime+1, 4) = MdrCasesActivated + MdrCasesEvolvedAct + MdrCasesEvolveAct_private; %active MDR disease incidence;
            diseaseCasesMat(timePeriod-burnInTime+1, 5) = sum(updating);
            diseaseCasesMat(timePeriod-burnInTime+1, 6) = MdrCasesTransmitted;
            diseaseCasesMat(timePeriod-burnInTime+1, 7) = hadTBbefore_ds;  %num ppl infected with DS that had TB before (reinfection case)
            diseaseCasesMat(timePeriod-burnInTime+1, 8) = hadTBbefore_mdr;  %num ppl infected with MDR that had TB before (reinfection case)
            
            %         deathsPostBurnIn(timePeriod-burnInTime+1, 1) = sum(deathsCounterMat(timePeriod,3:10)) - sum(deathsCounterMat(timePeriod,13:17)); %deaths from DS TB
            %         deathsPostBurnIn(timePeriod-burnInTime+1, 2) = sum(deathsCounterMat(timePeriod,13:17));  %deaths from MDR
            %         deathsPostBurnIn(timePeriod-burnInTime+1, 3) = numCurrentPpl;
            
            %get the dead ages
            runningDead = logical(runningDead);
            hasDSTB = logical(runningDSDead);
            hasMDRTB = logical(runningMDRDead);
            %%COMMENTED OUT FOR SPEED
            %         smoker = false(size(stateMat,1),1);
            %         %     hasDSTB(runningDead,1) = (stateMat(runningDead,Phealth) == 3 );
            %         %     hasMDRTB(runningDead,1) = (stateMat(runningDead,Phealth) == 4 );
            %         smoker(runningDead,1) = (stateMat(runningDead,Psmoking) == 1 );
            %         nosmokeDeathAgesMat(timePeriod-burnInTime+1,:) = histc(stateMat(runningDead & ~hasDSTB & ~hasMDRTB & ~smoker, Page), [0:1:100])';
            %         smokeDeathAgesMat(timePeriod-burnInTime+1,:)  = histc(stateMat(runningDead & ~hasDSTB & ~hasMDRTB & smoker, Page), [0:1:100])';
            %         nosmokeDSTBdeathAgesMat(timePeriod-burnInTime+1,:)  = histc(stateMat(runningDead & hasDSTB & ~smoker, Page), [0:1:100])';
            %         smokeDSTBdeathAgesMat(timePeriod-burnInTime+1,:)  = histc(stateMat(runningDead & hasDSTB & smoker, Page), [0:1:100])';
            %         nosmokeMDRTBdeathAgesMat(timePeriod-burnInTime+1,:)  = histc(stateMat(runningDead & hasMDRTB & ~smoker, Page), [0:1:100])';
            %         smokeMDRTBdeathAgesMat(timePeriod-burnInTime+1,:)  = histc(stateMat(runningDead & hasMDRTB & smoker, Page), [0:1:100])';
        end
    end
    %timer
    timeLoopTimer(timeLoopIndex) = timeLoopTimer(timeLoopIndex) + (cputime-timerStart);
    timeLoopTimerName{timeLoopIndex} = 'record';
    timeLoopIndex = timeLoopIndex + 1;
    timerStart = cputime;
    
    %     %%%%%%%%%%%%%%%%%%GETTING THE COUNTS FOR UPTAKE STUFF DELETE THIS DEBUGGING %%%%%%%%%
    %     if timePeriod == 1737
    % %         TBppl2010healthy = zeros(12,4);
    % %         TBppl2010sick = zeros(12,4);
    % %         TBppl2010diagnosed = zeros(12,4);
    % %         ageArray = [0, 20, 30, 40, 50, 60, 70, 101];
    % %         for sickness = 1:3
    % %             if sickness == 1
    % %                 sicknessVals = [3,4];
    % %                 sick = subsetterFunc(stateMat,[-1],Phealth, sicknessVals);
    % %             elseif sickness == 2
    % %                 sicknessVals = [0,1,2,3,4];
    % %                 sick = subsetterFunc(stateMat,[-1],Phealth, sicknessVals);
    % %             elseif sickness == 3
    % %                 sick = detectedPpl;
    % %             end
    % %
    % %             for agething = 1:7
    % %                 agelower = ageArray(agething);
    % %                 ageupper = ageArray(agething+1)-1;
    % %                 inAgePpl = subsetterFuncInequality(stateMat, sick,Page,agelower,ageupper);
    % %                 for demog = 1:4
    % %                     if demog == 1
    % %                         gender = 2;
    % %                         urbRur = 1;
    % %                     elseif demog == 2
    % %                         gender = 2;
    % %                         urbRur = 2;
    % %                     elseif demog == 3
    % %                         gender = 1;
    % %                         urbRur = 1;
    % %                     elseif demog == 4
    % %                         gender = 1;
    % %                         urbRur = 2;
    % %                     end
    % %                     sex = subsetterFunc(stateMat,inAgePpl,Psex,[gender]);
    % %                     urbRural = subsetterFunc(stateMat,sex,PurbanRur,urbRur);
    % %
    % %                     if sickness == 1
    % %                         TBppl2010sick(agething,demog) = size(urbRural,1);
    % %                     elseif sickness == 2
    % %                         TBppl2010healthy(agething,demog) = size(urbRural,1);
    % %                     elseif sickness == 3
    % %                         TBppl2010diagnosed(agething,demog) = size(urbRural,1);
    % %                     end
    % %                 end
    % %             end
    % %         end
    % %         TBppl2010healthy
    % %         TBppl2010sick
    % %         TBppl2010diagnosed
    % %
    %
    %         disp('total pop')
    %         numCurrentPpl
    %
    %         disp('diagnosisCounterMat')
    %         diagnosisCounterMat((timePeriod-3):timePeriod, :)
    %
    %         disp('treatment demog')
    %         treatmentDemog(1737-burnInTime ,:)
    %
    %     end
    
    %     %record age and health states
    %     if timePeriod >= burnInTime
    %         monthlyAgeHealthIndex = 1;
    %         for ageNum = 0:5:95
    %             for healthNum = 0:1:6
    %                 monthlyAgeHealth(timePeriod-burnInTime+1, monthlyAgeHealthIndex) = sum(stateMat(:,Phealth) == healthNum &  stateMat(:,Page) >= ageNum & stateMat(:,Page) < ageNum+5  );
    %                 monthlyAgeHealthIndex = monthlyAgeHealthIndex+1;
    %             end
    %         end
    %     end
    
    %extra validation measure made on Feb 14 2014
    if timePeriod == 1669
        disp('fast, slow, and old activation')
        TBparams.latentToActLessT2Yr
        TBparams.latentToActGreatT2Yr'
        latentGreaterT14Yrs_prob'
        totAct = TBparams.latentToActLessT2Yr(end, end) * ones(size(TBparams.latentToActLessT2Yr ,1), 200);
        totAct(1:size(TBparams.latentToActLessT2Yr ,1), 1:(size(TBparams.latentToActLessT2Yr ,2)+size( TBparams.latentToActGreatT2Yr,2)) ) = [TBparams.latentToActLessT2Yr, repmat(TBparams.latentToActGreatT2Yr, size(TBparams.latentToActLessT2Yr ,1), 1)];
        halfYearOldAct = reshape([latentGreaterT14Yrs_prob;latentGreaterT14Yrs_prob],2*size(latentGreaterT14Yrs_prob,2),1)';  %turn into half years
        
        for rowNum =1:size(TBparams.latentToActLessT2Yr)
            age = TBparams.latentToActAgeBrac(rowNum,1);
            endOfRow = min(age+60+size(halfYearOldAct,2)-1, size(totAct,2));
            totAct(rowNum,age+60:endOfRow) = halfYearOldAct + totAct(rowNum,age+60:endOfRow);
            if endOfRow ~= size(totAct,2)
                totAct(rowNum,endOfRow+1:end) = halfYearOldAct(end) + totAct(rowNum,endOfRow+1:end) ;
            end
        end
        tableHeader = 'time since infection 0,0.5, 1, 1.5, (rows are ages 0 5 10 etc)';
        tablePrinter(tableHeader, totAct, 'activationSurface', folderName);
        
        disp('transmission matrix')
        TBparams.modMonRateMat
    end
    
    if (timePeriod == 1669 || timePeriod == 1645 || timePeriod == 1585 || timePeriod == 1525) %Jan 2005, 2003, 1999, or 1993
        ageCounter = 0;
        for subsetAge = 0:5:95
            ageCounter = ageCounter + 1;
            inAge = subsetterFuncInequality(stateMat, find(updating==1),Page,subsetAge,subsetAge+4);
            privateAged = subsetterFuncInequality(stateMat, detectedPpl_private,Page,subsetAge,subsetAge+4);
            fastActAged = subsetterFuncInequality(stateMat, fastActivators,Page,subsetAge,subsetAge+4);%DEBUGGING
            slowActAged = subsetterFuncInequality(stateMat, slowActivators,Page,subsetAge,subsetAge+4);%DEBUGGING
            OLDslowActivatorsAged = subsetterFuncInequality(stateMat, OLDslowActivators,Page,subsetAge,subsetAge+4);%DEBUGGING
            
            for subsetSex = 1:2
                sexed = subsetterFunc(stateMat,inAge,Psex,[subsetSex]);
                hasTBin2003(ageCounter,subsetSex) = size(sexed,1);
                hasTB = subsetterFunc(stateMat,sexed,Phealth,[3,4]);
                hasTBin2003(ageCounter,2+subsetSex) = size(hasTB,1);
                hasTBtreatedB4 = subsetterFunc(stateMat,hasTB,PpastTreatment,[1,2,3]);
                hasTBin2003(ageCounter,4+subsetSex) = size(hasTBtreatedB4,1);
                hasTBinTreatment = subsetterFunc(stateMat,hasTB,Ptreatment,[1,2,4]);
                hasTBin2003(ageCounter,6+subsetSex) = size(hasTBinTreatment,1);
                privateSexed = subsetterFunc(stateMat,privateAged,Psex,subsetSex);
                hasTBin2003(ageCounter,8+subsetSex) = size(privateSexed,1);
                latTB = subsetterFunc(stateMat,sexed,Phealth,[1,2]);
                hasTBin2003(ageCounter,10+subsetSex) = size(latTB,1);
                
                fastActSexed = subsetterFunc(stateMat,fastActAged,Psex,subsetSex);      %DEBUGGING
                hasTBin2003(ageCounter,12+subsetSex) = size(fastActSexed,1);%DEBUGGING
                slowActSexed = subsetterFunc(stateMat,slowActAged,Psex,subsetSex);     %DEBUGGING
                hasTBin2003(ageCounter,14+subsetSex) = size(slowActSexed,1);%DEBUGGING
                OLDslowActivatorsSexed = subsetterFunc(stateMat,OLDslowActivatorsAged,Psex,subsetSex);     %DEBUGGING
                hasTBin2003(ageCounter,16+subsetSex) = size(OLDslowActivatorsSexed,1);%DEBUGGING
                
            end
        end
        %print out
        tableHeader = 'male, female, male hasTB, female hasTB, male hasTBinTrt, female hasTBinTrt, male hasTBhadTrt, female hasTBhadTrt, (rows are ages 0 5 10 etc)';
        tablePrinter(tableHeader, hasTBin2003, strcat('hasTBin_',num2str(timePeriod)), folderName);
    end
    
    timePeriod = timePeriod + 1;  %increment time period, since now it's a while loop and needs incrementation.
    if LEbuilder == 1
        if timePeriod == totPeriods + 1  %reached the end of time for this LEbuilder cohort, plus 1 since just incremented
            
            %%%%%%%%%%%% write down what happened to this cohort %%%%%%%%%%%%
            %age, sex, health, treatment status
            lastState = zeros(1,20);
            counter = 1;
            for i = [0,1,2,3,4]
                healthed = subsetterFunc(stateMat,[-1],Phealth,[i]);
                for j = [0,1,2,4]
                    lastState(counter) = size(subsetterFunc(stateMat,healthed,Ptreatment,j),1);
                    counter = counter+1;
                end
            end
            
            %reset time
            disp(cohortNum)
            timePeriod = origTotPeriods+1;  %reset the time periods for next cohort
            if cohortNum == totLEbuilderCohorts  %unless just finished the last cohort, in which case really finished
                tablePrinter('undiscounted costs incl burninTime', LEcostsMatBig, strcat('costsOverTime_allCohorts'), folderName);
                tablePrinter('over health 012346 and trt 0124', lastStateBig, strcat('lastState_allCohorts'), folderName);
                tableHeader = 'healthyLatNeverTrt,healthyLatPastTrt,healthyLatDOTS,healthyLatMDRTrt,DSpplNoTrt,DSpplInDOTS,DSpplInMDR,MDRpplNotInMDR,MDRpplInMDR';
                tablePrinter(tableHeader, LEnumQalyPpl_males', 'numQalyPpl_males', folderName);
                tablePrinter(tableHeader, LEnumQalyPpl_females', 'numQalyPpl_females', folderName);
                timePeriod = totPeriods+1;
            end
            %print LE costs and QALYs
            LEcostsMatBig(:,cohortNum) = LEcostsMat;
            lastStateBig(cohortNum,:) = lastState;
            numMaleCohorts = totLEbuilderCohorts/2;
            if cohortNum <= numMaleCohorts
                LEnumQalyPpl_males(:,(cohortNum*LEbuilderCohortDuration)-(LEbuilderCohortDuration-1):(cohortNum*LEbuilderCohortDuration)) = LEnumQalyPpl; %for LEbuilder
            else
                LEnumQalyPpl_females(:,((cohortNum-numMaleCohorts)*LEbuilderCohortDuration)-(LEbuilderCohortDuration-1):((cohortNum-numMaleCohorts)*LEbuilderCohortDuration)) = LEnumQalyPpl; %for LEbuilder
            end
            
        end
    end
    
end  %%end time loop
toc
fprintf('Loop timer statistics:\n');
fprintf('   %20s | %7s\n','Description','Time (s)');
fprintf('   ---------------------+---------\n');
for ii = 1:length(timeLoopTimer)
    fprintf('   %20s | %7.3f\n',timeLoopTimerName{ii},timeLoopTimer(ii));
end

%%%%%%%%%%%% write down end stateMat %%%%%%%%%%%%
%age, sex, health, treatment status

lastState = zeros(1,4000);
counter = 1;
for sex = 1:2
    sexed = subsetterFunc(stateMat,[-1],Psex,[sex]);
    for i = 0:4
        healthed = subsetterFunc(stateMat,sexed,Phealth,[i]);
        for j = [0,1,2,4]
            trted = subsetterFunc(stateMat,healthed,Ptreatment,[j]);
            for age = 0:1:99
                lastState(counter) = size(subsetterFunc(stateMat,trted,Page,[age]),1);
                %lastStateCat(counter) = strcat(num2str(sex),num2str(i),num2str(j),num2str(age));
                counter = counter+1;
            end
        end
    end
end

%%%%%%%%%%%% DECOMPRESSIONS %%%%%%%%%%%%
if (recordUnit == 1 && totPeriods >= burnInTime && LEbuilder == 0)
    
    disp('decompressing reincarnations')
    disp('size of recNo is')
    size(recNo)
    
    %make matrix with only the Phealth values
    healthMat = floor( (  mod(recNo,48) )/6 );
    
    % Decompress the reincarnation thing
    [reincRows,reincCols] = find(healthMat == 6);  %the 6's are all dead
    % Sort them by row, then by columns in case of tie.
    % We will implement this by sorting by column first,
    % and then by row (each sort preserves the old order in case of tie)
    [reincCols,sortIndex] = sort(reincCols);
    reincRows = reincRows(sortIndex);
    %% sort rows and cols based on rows
    [reincRows,sortIndex] = sort(reincRows);
    reincCols = reincCols(sortIndex);
    
    numReincarnations = size(reincRows,1);
    % Start with the states from before, using the compressed recNo matrix
    decompressedStates = zeros(    size(recNo,1)+ numReincarnations  ,    size(recNo,2)    );
    decompressedStates(1:size(recNo,1),:) = recNo;
    
    disp('building decompressedStates')
    fprintf('total reborn index nums %i\n',numReincarnations);
    % Copy the reincarnations to the bottom of this array
    for rebornIndex = 1:numReincarnations
        if(mod(rebornIndex,10000)==0 )
            fprintf('Doing reborn index %i out of %i\n',rebornIndex,numReincarnations);
        end
        % Read in the next reincarnation event
        originalIndex = reincRows(rebornIndex);
        reincarnationTime = reincCols(rebornIndex);
        % Figure out the end of this life history
        if (rebornIndex < numReincarnations && reincRows(rebornIndex+1) == originalIndex)
            % If this person got reincarnated again, only copy up to the next
            % reincarnation
            stopTime = reincCols(rebornIndex+1);
        else
            % Else, this person lived up to the end of the simulation
            stopTime = size(recNo,2);
        end
        % Fill the beginning of the new life history with not yet born
        decompressedStates(size(stateMat,1)+rebornIndex,1:reincarnationTime) = 87342;  %stateMat translation: [1,0,0,0,2,0,5,0,0,0,0,0]
        % Copy over the life history
        decompressedStates(size(stateMat,1)+rebornIndex,reincarnationTime+1:stopTime) =....
            decompressedStates(originalIndex,reincarnationTime+1:stopTime);
        % Erase the life history from the old location (fill with dead)
        decompressedStates(originalIndex,reincarnationTime+1:stopTime) = ...
            decompressedStates(originalIndex,reincarnationTime);
        % Fill the end of the newly copied life history (if they died again) with
        % dead
        decompressedStates(size(stateMat,1)+rebornIndex,stopTime+1:end) = ...
            decompressedStates(size(stateMat,1)+rebornIndex,stopTime);
    end
    %rewrite states
    recNo = decompressedStates;
    
    % bugcatcher code
    for timePeriod = 1:size(recNo,2)
        sMatDec = recNoDecoder(recNo(:,timePeriod));
        if any(sMatDec(:,Psex) > 2 | sMatDec(:,Psex) < 1)
            warning('Decoder found sex out-of-bounds');
        end
        if any(sMatDec(:,Page) > 99 | sMatDec(:,Page) < 0)
            warning('Decoder found age OOB');
        end
        if any(sMatDec(:,Psmoking) > 1 | sMatDec(:,Psmoking) < 0)
            warning('Dec found smoking OOB');
        end
        if any(sMatDec(:,PurbanRur) > 2 | sMatDec(:,PurbanRur) < 1)
            warning('Dec found urban/rural OOB');
        end
        if any(sMatDec(:,Ptreatment) > 4 | sMatDec(:,Ptreatment) < 0 ...
                | sMatDec(:,Ptreatment) == 3)
            warning('Dec found treatment OOB');
        end
        if any(sMatDec(:,Phealth) > 6 | sMatDec(:,Phealth) < 0 )
            warning('Dec found health OOB');
        end
    end
    
end
if (totPeriods > burnInTime && recordUnit == 1)
    %intRecNoPostBurnIn = int32(recNo(:,((burnInTime/recordUnit)+1):(totPeriods/recordUnit)));
    intRecNo = int32(recNo);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%make graphs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if writeBurnIn ~= 1
    %  for natural death rates debug and calibration
    if durationYrs == 230
        currentDirectory = pwd;
        cd(folderName);
        MonSinceActTruncations = sparse(MonSinceActTruncations);
        MonSinceAct = sparse(MonSinceAct);
        MonSinceLat = sparse(MonSinceLat);
        save fishMatrices.mat deathAgesMat aliveAgesMat aliveActTBAgesMat aliveActTBAgesMat_after MonSinceActTruncations MonSinceAct MonSinceLat
        cd(currentDirectory);
    end
    
    %force of infection
    tableHeader = 'effective monthly contact rate, recorded once every year';
    tablePrinter(tableHeader, forceOfInfection, 'forceOfInfection', folderName);
    
end


if (strcmp(plotResolution,'-r0') == 0 && LEbuilder == 0)
    
    %state at the end of the simulation
    tablePrinter('see code lastState', lastState, 'lastState', folderName);
    
    %QALY tables
    %tableHeader = 'healthy_noTreat, healthy_noTreat_pastTreat, latDS_noTreat, latDS_noTreat_pastTreat, latMDR_noTreat, latMDR_noTreat_pastTreat, actDS_noTreat, actDS_catI, actDS_catII, actDS_catIII, actDS_catIV, actMDR_noTreat, actMDR_catI, actMDR_catII, actMDR_catIII, actMDR_catIV, notYetBorn, dead';
    tableHeader = 'healthyLatNeverTrt,healthyLatPastTrt,healthyLatDOTS,healthyLatMDRTrt,DSpplNoTrt,DSpplInDOTS,DSpplInMDR,MDRpplNotInMDR,MDRpplInMDR';
    tablePrinter(tableHeader, numQalyPpl', 'numQalyPpl_postBurnIn', folderName);
    %costs tables
    tableHeader = 'undiscounted costs incl burninTime';
    tablePrinter(tableHeader, costsMat, 'costsOverTime', folderName);
    
    %     if totPeriods >= burnInTime
    %         %QALYS csv
    %         tablePrinter('totalQalys_byMonth', Qalys, 'TotalQALYs_scaledToIndianPop', folderName);
    %
    %         %deaths csv
    %         tablePrinter('numDSdeaths, numMDRdeaths', popDeaths, 'TotalDeaths_scaledToIndianPop', folderName);
    %
    %         %cases csv
    %         tablePrinter('latent DS infection incidence,latent MDR infection incidence,active DS disease incidence,active MDR disease incidence', popCases, 'TotalCases_scaledToIndianPop', folderName);
    %     end
    
    
    %monthly actTB table
    tableHeader = 'healthy, latDS, latMDR, DSnoTreat, DScatI, DScatII, DScatIV, MDRnoTreat, MDRcatI, MDRcatII, MDRcatIV, numCurrentPpl, only accurate if RNTCP spec is 1';
    tablePrinter(tableHeader, monthlyActOutcomes, 'monthlyActOutcomes', folderName);

    %numPplBorn
    tableHeader = 'numPpl born in each time period';
    tablePrinter(tableHeader, numPplBorn, 'numPplBorn', folderName);

    if TBparams.limitOutputs ~= 1
        %deaths graph
        %    csvDeathsPlotter(deathsCounterMat, folderName, plotResolution);
        tableHeader = 'nonsmoking natural deaths, smoking natural deaths, nonsmoking untreatedTB deaths, smoking untreatedTB deaths, catI nonsmoking deaths, catI smoking deaths, catII nonsmoking deaths, catII smoking deaths, catIV nonsmoking deaths, catIV smoking deaths, nonsmoking age-related deaths (reach age 99), smoking age-related deaths (reach age 99), all active mdrTB related deaths (see code), all active sensTB related deaths, col 19 is num private deaths, col23 is TBppl who died from nonTB causes';
        tablePrinter(tableHeader, deathsCounterMat, 'deathsCounterMat', folderName);
        
        %counting TB deaths and self cure
        %     tableHeader = 'numUntrtedTBdeaths, numTBselfCure';
        %     tablePrinter(tableHeader, trtCountMat, 'trtCountMat', folderName);
        
        %deaths to make the deaths averted numbers for different treatments
        tableHeader = 'numDSdeaths, numMDRdeaths, numCurrentPpl';
        tablePrinter(tableHeader, deathsPostBurnIn, 'TBdeaths_postBurnIn', folderName);
        
        %deathAges tables
        tableHeader = 'num dead in each age (time period by rows and each column is an age 0 to 100 inclusive)';
        tablePrinter(tableHeader, nosmokeDeathAgesMat, 'nosmokeDeathAgesMat', folderName);
        tablePrinter(tableHeader, smokeDeathAgesMat, 'smokeDeathAgesMat', folderName);
        tablePrinter(tableHeader, nosmokeDSTBdeathAgesMat, 'nosmokeDSTBdeathAgesMat', folderName);
        tablePrinter(tableHeader, smokeDSTBdeathAgesMat, 'smokeDSTBdeathAgesMat', folderName);
        tablePrinter(tableHeader, nosmokeMDRTBdeathAgesMat, 'nosmokeMDRTBdeathAgesMat', folderName);
        tablePrinter(tableHeader, smokeMDRTBdeathAgesMat, 'smokeMDRTBdeathAgesMat', folderName);
        tableHeader = 'MDR_cat1Dead,DS_cat1Dead,totalCat1,MDRinCat1';
        tablePrinter(tableHeader, cat1MDRdead, 'cat1MDRdead', folderName);
        
        tableHeader = 'MDR_cat2Dead,DS_cat2Dead,totalCat2,MDRinCat2';
        tablePrinter(tableHeader, cat2MDRdead, 'cat2MDRdead', folderName);
        
        tableHeader = 'total to private, num stay private, num to RNTCP';
        tablePrinter(tableHeader, numRefPpl, 'PPMrefPpl', folderName);
        
        tableHeader = 'total death ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, deathAgesMat_TBmac , 'deathAgesMat_TBmac', folderName);
        
        tableHeader = 'TB death ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, deathTBAgesMat_TBmac , 'deathTBAgesMat_TBmac', folderName);

        tableHeader = 'Alive ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, aliveAgesMat_TBmac , 'aliveAgesMat_TBmac', folderName);
        
        tableHeader = 'active TB Alive ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, aliveTBAgesMat_TBmac , 'aliveTBAgesMat_TBmac', folderName);
        
        tableHeader = 'active MDR TB Alive ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, aliveMDRAgesMat_TBmac , 'aliveMDRAgesMat_TBmac', folderName);
        
        tableHeader = 'lat TB Alive ages over time post burn in.  Cols are ages 0-9 10-15 16-19 20-29 30-39 etc to 90-109';
        tablePrinter(tableHeader, aliveLatAgesMat_TBmac , 'aliveLatAgesMat_TBmac', folderName);
        
        tableHeader = 'Number of diagnosed in cat I II not double counting, mdr tested, mdr tested pos, diagnosedPrivate, adults diagnosed catI and II, adults diagnosedPrivate, trtNaive private diagnoses (not notifications)';
        tablePrinter(tableHeader, notification_TBmac , 'notification_TBmac', folderName);
        
        tableHeader = 'numScreened intervention ACF, num diagnosed by intervention ACF';
        tablePrinter(tableHeader, numScreened , 'numScreened_TBmac', folderName);
        
        %DOTS performance measures table
        tablePrinter('numCATI, numCATIdefaults, numCATIgraduating, numCATIcured, numCATII, numCATIIdefaults, numCATIIgraduating, numCATIIcured, numCATIV, CATIVgraduating, CATIVcured', DOTSperformanceMeas, 'DOTSperformanceMeas', folderName);
        
        %past treatment counter
        tablePrinter('actTBPpastTreatment_0 , actTBPpastTreatment_1, actTBPpastTreatment_2, actTBPpastTreatment_3, actTBAnyPastPrivate', pastTreatment, 'pastTreatmentCounts', folderName);
        %catI types of nonsuccess
        tablePrinter('defaults,  success_but_latent, failures', catInonSuccessTypes , 'catInonSuccessTypes', folderName);
        
        %catII types of nonsuccess
        tablePrinter('eightMonthers, DS_8monthers, gotLatMDR, MDRcured, DS_failed, failedGotMDR', catIInonSuccessTypes , 'catIInonSuccessTypes', folderName);
        
        %treatment demogs table
        tablePrinter('catI_male, catI_youngMale, catI_female, catI_youngFemale, catII_male, catII_youngMale, catII_female, catII_youngFemale, numCurrentPpl', treatmentDemog, 'treatmentDemog', folderName);
        
        %average age of smokers and nonsmokers
        if durationYrs == 230
            tablePrinter('smokers ave age, nonsmokers ave age, smokers Post30 aveAge, nonsmokers Post30 aveAge', aveSmokerNonAge, 'aveAge_smokersNonsmokers', folderName);
        end
        
        %health states by treatment category tables
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', NoTreatHealthStates, 'healthByTreat_NoTreatHealthStates', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIHealthStates, 'healthByTreat_CatIHealthStates', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIIHealthStates, 'healthByTreat_CatIIHealthStates', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIVHealthStates, 'healthByTreat_CatIVHealthStates', folderName);
        %
        %     %updating health states by treatment category tables
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', NoTreatHealthStates_updating, 'healthByTreat_NoTreatHealthStates_updating', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIHealthStates_updating, 'healthByTreat_CatIHealthStates_updating', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIIHealthStates_updating, 'healthByTreat_CatIIHealthStates_updating', folderName);
        %     tablePrinter('health0, health1, health2, health3, health4, health5, health6, extra', CatIVHealthStates_updating, 'healthByTreat_CatIVHealthStates_updating', folderName);
        
        
        tableHeader = 'numPpl_1privateClinic, numPpl_2privateClinic, numPpl_3privateClinic,etc';
        tablePrinter(tableHeader, numPrivateClinics, 'numPrivateClinics', folderName);

        tableHeader = 'average number of months before got on treatment, died, or self cured';
        tablePrinter(tableHeader, ignoredTBmat, 'aveMonIgnoredTBmat', folderName);

        tableHeader = 'num MDR activations post burnIn';
        tablePrinter(tableHeader, MdrActivation, 'numMDRactivations_postBurnIn', folderName);
        
        tableHeader = 'num cases MDR evolved, num cases MDR transmitted';
        tablePrinter(tableHeader, MdrMethod, 'numNewMDR_evolvedOrTrans_postBurnIn', folderName);
        
        %     tableHeader = 'numPpl age0to4 health0, age0to4 healthstate1,age0to4 health2,age0to4 healthstate3,age0to4 health4,age0to4 health5,age0to4 health6, age5to9 health0, age5to9 health1,';
        %     tablePrinter(tableHeader, monthlyAgeHealth, 'monthlyAgeHealth', folderName);
        
        
        tableHeader = 'num kids 1 to 9 healthy, latDS, latMDR, actDS, actMDR';
        tablePrinter(tableHeader, kids_TBstateCount, 'kids_TBstateCount', folderName);
        
        %TB cases for comparing across treatments for cases averted
        tableHeader = 'latent DS infection incidence,latent MDR infection incidence,active DS disease incidence,active MDR disease incidence,numCurrentPpl, num latent MDR cases transmitted (subset of col 2), number of DS infections had TB before (reinfection case), number of MDR infections had TB before (reinfection case)';
        tablePrinter(tableHeader, diseaseCasesMat, 'DSmdrIncidence_postBurnIn', folderName);
        
        %activation Table
        tableHeader = 'totalLatentLessT2yrs, totalLatentGreaterT2yrs, fastNonsmokerActivations, fastSmokerActivations, slowNonsmokerActivations, slowSmokerActivations, numunder15Activated, 15andOverFastAct';
        tablePrinter(tableHeader, activaCounterMat, 'latentToActive', folderName);
        tableHeader = 'male nonSmoker DS fast activations,male nonSmoker DS slow,male smoker DS fast,male smoker DS slow,female nonSmoker DS fast,female nonSmoker DS slow,female smoker DS fast,female smoker DS slow,male nonSmoker MDR fast,male nonSmoker MDR slow,male smoker MDR fast,male smoker MDR slow,female nonSmoker MDR fast,female nonSmoker MDR slow,female smoker MDR fast,female smoker MDR slow';
        tablePrinter(tableHeader, activaCountSex, 'latentToActive_sex', folderName);
        tableHeader = 'male nonSmoker DS fast latents,male nonSmoker DS slow latents,male smoker DS fast,male smoker DS slow,female nonSmoker DS fast,female nonSmoker DS slow,female smoker DS fast,female smoker DS slow,male nonSmoker MDR fast,male nonSmoker MDR slow,male smoker MDR fast,male smoker MDR slow,female nonSmoker MDR fast,female nonSmoker MDR slow,female smoker MDR fast,female smoker MDR slow';
        tablePrinter(tableHeader, latentCount, 'latentToActiveEligible_sex', folderName);
        
        %diagnosis tables
        tableHeader = 'totalPop,totActTB,actTBnoTreatment,diagnosedDOTS,diagnosedMDR, diagnosedToCatI, diagnosedToCatII, aveTimeSinceActivation, MDRdiagnosedToCatI, MDRdiagnosedToCatII,relapseDiagnosedToCatII,MDRrelapseDiagnosedToCatII';
        tablePrinter(tableHeader, diagnosisCounterMat, 'diagnosedPpl', folderName);
        %     tableHeader = 'age at diagnosis 0-4, 5-9, 10 - 14, etc';
        %     tablePrinter(tableHeader, diagnosedAges, 'diagnosedAges', folderName);
        tableHeader = 'timeFromAct_all, median_timeFromAct, ave num mons since act(men age 0-9), median(men age 0-9), num in grp(men age 0-9), ave num mons since act(females age 0-9), median(females age 0-9), num in grp(females age 0-9), same pattern age 10-19, etc';
        tablePrinter(tableHeader, diagnosisTimeMat, 'diagnosisTime_CatI', folderName);
        %     tableHeader = 'numPpl diagnosed with activation 1 month ago, 2months ago, 3months ago...';
        %     tablePrinter(tableHeader, diagnosisTimeDist, 'diagnosisCATI_Distribution', folderName);
        tableHeader = 'activeAll Ppl, coveraledAllPpl_maySeek, coveredAllPpl, canBeDetected, detectedByTest, actuallyDetected';
        tablePrinter(tableHeader, diagnosedPool, 'diagnosedPool', folderName);
        %
        %     tableHeader = 'healthyDSppl_MDRtest,healthyDSppl_MDRtestPos,mdrPpl_MDRtest,mdrPpl_MDRtestPos, MDRspecPool_healthy, MDRspecPool_healthyNoTrt, MDRspecPool_healthyDOTS, aveTimeinTrt, MDRspecPool_DS, MDRspecPool_DSNoTrt, MDRspecPool_DSDOTS, aveTimeinTrt,firstMonthDST_healthy, firstMonthDST_DS, firstMonthDST_MDR, catI_healthy, catI_DS, catI_mdr, catIImid_healthy, catIImid_DS, catIImid_mdr, catIIfin_healthy, catIIfin_DS, catIIfin_mdr';
        %     tablePrinter(tableHeader, mdrTestPpl, 'mdrTestPpl', folderName);
        
        
        %private treat
        tableHeader = 'totPeople, total actTBppl, tot actTB not on Trt, num ppl diagnosed private, ace mon since act';
        tablePrinter(tableHeader, privateCounterMat, 'privateCounterMat', folderName);
        
        %age graphs
        disp('AgeStructure')
        graphTable = plotByTime2(recNo, 0, 0, plotUnit,1);
        print('-dpng',plotResolution,[folderName '/AgeStructure.png']);
        tablePrinter(graphTable{1}, graphTable{2}, 'AgeStructure', folderName);
        close all
        
        disp('ActiveTB_AgeStructure')
        graphTable = plotByTime3(recNo, 0, 0, plotUnit,1);
        print('-dpng',plotResolution,[folderName '/ActiveTB_AgeStructure.png']);
        tablePrinter(graphTable{1}, graphTable{2}, 'ActiveTB_AgeStructure', folderName);
        close all
        
        %snapshot of ages at last time period
        disp('making outcome by age plot')
        stateMatR = recNoDecoder(recNo(:,(floor(totPeriods/recordUnit))));
        inithealthMat2(:,numHealthStates+2) = stateMatR(:,Page);
        for i = 0:numHealthStates
            inithealthMat2(:,i+1) = (stateMatR(:,Phealth) == i);
        end
        plotByAge(inithealthMat2, numHealthStates+2, [1 2 3 4 5 6 7]);
        title('Health Outcomes By Age');
        print('-dpng',plotResolution,[folderName '/healthOutcomeByAge.png']);
        close all
        
        
        %make the non-plotByTime graphs
        csvPlotMatrix(recNo, folderName);  %make 'csvMatForPlotting.txt';
        csvPlotter('csvMatForPlotting.txt', folderName, plotResolution);
        
        %make the two-way comparison graphs for the results
        %    resultsPlotter(recordUnit, folderName, plotResolution)
        
    end
    %incidence graph
    %incidenceMat(1:2,1:2) = [0,0;0,0];   %set the first two rows (time periods) of incidence to zero, so the scaling will not be too big
    csvIncidencePlotter(incidenceMat, folderName, plotResolution);
    tableHeader = 'activations, totalPop';
    tablePrinter(tableHeader, incidenceMat, 'incidenceMatrix_forMakingWHOComparison', folderName);
    
    
end
fprintf('Total time to run: %0.3f sec\n',cputime-startTime);
%========================
%close figures
close all
%Write to diary
if (recordUnit == 1 && totPeriods > burnInTime && durationYrs == 230)
    disp('writing to diary')
    filename = 'recNo.mat';
    save( [folderName '/' filename], 'intRecNo','-v7.3');  %intRecNoPostBurnIn
end
diary off
if (recordUnit == 1 && durationYrs == 230)
    clearvars -except folderName recNo
    
    %disp('running latent to active analysis since record unit is 1 ');
    %latToActAnalysis(folderName);
    
    %disp('building stata output');
    %stataOutputMaker(recNo, folderName);
end
disp('simulation complete')
