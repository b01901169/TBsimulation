%make discounted costs and qalys
%script used in qalyMaker2

%make discount factor


%scale up to Indian population levels and discount
%numPplHealthState =( simulatedHealthStates(startTime:endTime,1:end-1) ./ repmat(numSimPpl,1,(size(simulatedHealthStates,2)-1)) ) .* repmat(correctPopScale,1,(size(simulatedHealthStates,2)-1));
perPersonUtility =( simulatedHealthStates(startTime:endTime,1:end-1) * utilityWeights) ./ numSimPpl;
nondiscountedQALYs = sum(perPersonUtility);
discountedQALYs = correctDiscountFactor*(perPersonUtility);

%costs
aveCostPerPerson = simulatedCosts(startTime:endTime)./numSimPpl;
%scaledCosts = aveCostPerPerson .* correctPopScale;
nondiscountedCosts = sum(aveCostPerPerson);
discountedCosts = correctDiscountFactor*aveCostPerPerson;


