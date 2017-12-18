function [discountedCosts,discountedQALYs,cohortAges, LEbuilderNumCohorts] = qalyMaker2(start, endT, subFolder, LEbuilder,discountIndex, LEbuilderN)
%qalyMaker script makes qalys over the specified duration discounted to
%start time.  LEbuilderN is the number of people in the LEbuilder
%simulation run.

cohortAges=0;
LEbuilderNumCohorts = 0;

 dir = pwd;
% fileFolder = strcat(dir,'\',subFolder);
fileFolder = subFolder;

if (start >= 1996 && start <= 2126)
    startTime = ((start-1996)*12)+1;
else
    disp('Script does not handle this start time')
end

if (endT >= 1996 && endT <= 2126)
    endTime = ((endT-1996)*12)+1;
else
    disp('Script does not handle this end time')
end
   
if LEbuilder == 1 && (discountIndex >= 1996 && discountIndex <= 2126)
    discountI = ((discountIndex-1996)*12)+1;
elseif (LEbuilder == 1)
    disp('Script does not handle this discountIndex time')
end

%healthyLatNeverTrt	healthyLatPastTrt	healthyLatDOTS	healthyLatMDRTrt	DSpplNoTrt	DSpplInDOTS	DSpplInMDR	MDRpplNotInMDR	MDRpplInMDR
utilityWeightsYr = [...
1
0.941699
0.851699
0.806699
0.663
0.843
0.753
0.663
0.753
];

utilityWeights = utilityWeightsYr/12;  %since monthly

IndianPop = [...
1950	371857000	;...
1955	406374000	;...
1960	447844000	;...
1965	496400000	;...
1970	553874000	;...
1975	622097000	;...
1980	700059000	;...
1985	784491000	;...
1990	873785000	;...
1995	964486000	;...
2000	1053898000	;...
2005	1140043000	;...
2010	1224614000	;...
2015	1308221000	;...
2020	1386909000	;...
2025	1458958000	;...
2030	1523482000	;...
2035	1579802000	;...
2040	1627029000	;...
2045	1664519000	;...
2050	1692008000	;...
];

%Load historic population in India
%TBparams.IndianPop2 = interp1(IndianPop(:,1), IndianPop(:,2), [1996:1/12:2046])';
TBparams.IndianPop2 = interp1(IndianPop(:,1), IndianPop(:,2), [1996:1/12:2050])';
TBparams.IndianPop2 = [TBparams.IndianPop2; IndianPop(end,2)*ones(12*(2126-2050+1),1)]; %just use a constant population after 2050

%make discount rate 
TBparams.discountRate = 0.03;
%TBparams.discountRate = 0;
discountFactor = (1/(1+TBparams.discountRate)).^([0:1/12:(130)]);
% expVec = reshape(repmat([0:1:130],12,1) .* ones(12,131),12*131,1);
% discountFactor = (1/(1+TBparams.discountRate)).^(expVec);


%get health states from simulation
cd(fileFolder);
if LEbuilder == 0
    simulatedHealthStates = dlmread('numQalyPpl_postBurnIn.csv', ',' ,1,0);
    simulatedHealthStates = [simulatedHealthStates(:,1:9),simulatedHealthStates(:,19)];  %first 9 columns and the 19th col are used
    simulatedCosts = dlmread('costsOverTime.csv', ',' ,1,0);
    numSimPpl = simulatedHealthStates(startTime:endTime,end);
    correctPopScale = TBparams.IndianPop2(startTime:endTime,1);
    correctDiscountFactor = discountFactor(1,1:endTime-startTime+1);    
    
    cd(dir)
    costAndQalyMaker;
    %   discountedQALYs,discountedCosts
        
else  %doing LEbuilder files
    cd(fileFolder);
    simulatedCosts_allCohorts = dlmread('costsOverTime_allCohorts.csv', ',' ,1,0);
    simulatedHealthStates_males = dlmread('numQalyPpl_males.csv', ',' ,1,0);
    simulatedHealthStates_females = dlmread('numQalyPpl_females.csv', ',' ,1,0);
    simulatedHealthStates_all = [simulatedHealthStates_males;simulatedHealthStates_females];

    cd(dir);  %now need to change to outputs folder, which should be one up
    cd ..      
    cohortValues = dlmread('LEbuilderStates_values.csv', ',' ,1,0);
    index = 1; sex = 2; healthState= 3; trtState = 4; age = 5; %columns
    cohortValues = [[1:1:size(cohortValues,1)]',cohortValues];  %make a matrix of the values (index, sex, health, trt, age)
    cohortAges = unique(sort(cohortValues(:,age),'descend') );
    maxCohortAge = max(cohortAges);
    ageIncr = [cohortAges(2:end);100]-cohortAges;

%     disp('LEbuilderCohortDuration')
    LEbuilderCohortDuration = size(simulatedCosts_allCohorts,1);
    LEbuilderNumCohorts = size(simulatedCosts_allCohorts,2);
    LEqalyCostHorizonTable = zeros(LEbuilderNumCohorts,3);
    numSimPpl = LEbuilderN*ones(LEbuilderCohortDuration,1); %only one year at a time now.  %numSimPpl = LEbuilderN*ones(endTime-startTime,1);  %no births, always 10000 people scale factor.  changed to LEbuilderN number of ppl.
    %correctPopScale = TBparams.IndianPop2(startTime:endTime-1,1);
    yearsLived =[1/24:1/12:(LEbuilderCohortDuration/12)-1/24]';
    
    for LEcohortNum = LEbuilderNumCohorts:-1:1 %1:400 400 cohorts.  actually more cohorts now since doing all ages
    %for LEcohortNum = 3:-1:1 % debugging
        simulatedCosts = simulatedCosts_allCohorts(:,LEcohortNum);
        simHealthval = LEbuilderCohortDuration*LEcohortNum;
        simulatedHealthStates = simulatedHealthStates_all(simHealthval-(LEbuilderCohortDuration-1):simHealthval,[1:9, 19]);%first 9 columns and the 19th col are used 
      
        %make conditional life expectancy       
        numAlive = sum(simulatedHealthStates(:,1:end-1),2);
        numAlivePriorPeriod = [LEbuilderN;numAlive(1:end-1)];
        numDied = numAlivePriorPeriod-numAlive;
        totSurvivedWholePeriod = LEbuilderN - sum(numDied);
        totYrsLived = (totSurvivedWholePeriod*(LEbuilderCohortDuration/12)) + (numDied(1:end)' * yearsLived(1:end));
        condLEqalyCostHorizonTable(LEcohortNum,1) = totYrsLived/LEbuilderN;  %first col is LE
        
        %make cost and qalys discounted to the first year that cohort was simulated for
        cd(dir)
        endTime = size(simulatedHealthStates,1);
        %costAndQalyMaker;  %not need the same kind of discounting now, just do it here
        perPersonUtility =( simulatedHealthStates(1:endTime,1:end-1) * utilityWeights) ./ numSimPpl;
        partiallyDiscountedQALYs = discountFactor(1:LEbuilderCohortDuration)*perPersonUtility;% sum(perPersonUtility);
        condLEqalyCostHorizonTable(LEcohortNum,3) = partiallyDiscountedQALYs; % discountedQALYs;  %3nd col is qalys
        
        %costs
        aveCostPerPerson = simulatedCosts(1:endTime)./numSimPpl;
        partiallyDiscountedCosts = discountFactor(1:LEbuilderCohortDuration)*aveCostPerPerson;% sum(aveCostPerPerson);
        condLEqalyCostHorizonTable(LEcohortNum,2) = partiallyDiscountedCosts; % discountedCosts;  %2rd col is costs
        
%         if mod(LEcohortNum,100) == 0
%             disp(num2str(LEcohortNum));
%         end
    end
    
    %collapse: start at end and roll up
    cd(fileFolder);
    weightMat = dlmread('lastState_allCohorts.csv', ',' ,1,0);
    weightMat = weightMat/LEbuilderN; %normalize weights to prob mat
    
    for sexVal = 1:2
        for outputTypes = 1:3  %for LE, costs, and QALYs
            correctRows99yos = find(cohortValues(:,sex) == sexVal & cohortValues(:,age) == maxCohortAge  );  %initialize with 99 year olds.  now using the oldest cohort may not be 99
            LEqalyCostHorizonTable(correctRows99yos,outputTypes) = condLEqalyCostHorizonTable(correctRows99yos,outputTypes);
            
            for ageind = ( size(cohortAges,1)-1 ):-1:1
                ageVal = cohortAges(ageind);                
                correctRowsNextYear = find(cohortValues(:,sex) == sexVal & cohortValues(:,age) == ageVal+ageIncr(ageind) );
                correctRowsThisYear = find(cohortValues(:,sex) == sexVal & cohortValues(:,age) == ageVal );
                % cohortValues(correctRowsNextYear,:);  %can check that these are the health / trt states that you want
                discountFac = (1/(1+TBparams.discountRate))^(LEbuilderCohortDuration/12);                
                if outputTypes == 1
                    discountFac = 1;
                end
                LEqalyCostHorizonTable(correctRowsThisYear,outputTypes) = (discountFac*weightMat(correctRowsThisYear,1:20)*LEqalyCostHorizonTable(correctRowsNextYear,outputTypes)) + condLEqalyCostHorizonTable(correctRowsThisYear,outputTypes);  %conditional life years
            end
        end
    end
    %discount costs and QALYs to the correct start time (discounted to the first year of first cohort at this point)    
    LEqalyCostHorizonTable(:,2:3) = ( (1/(1+TBparams.discountRate))^10 )* LEqalyCostHorizonTable(:,2:3);  %discount by the 10-year-base analysis period

    %put the table in the discountedCosts output and make discounted Qalys empty to use the same names as nonLEbuilder run
    discountedCosts = LEqalyCostHorizonTable;
    discountedQALYs = [];
    
end
cd(dir);





