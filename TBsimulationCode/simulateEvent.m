function outcome = simulateEvent(prevalence, default, ageBracTable, age, sex)
% outcome = simulateEvent(prevalence, default, ageBrac, age, sex)
%
% Function to simulate an outcome based on a prevalence table
% Returns:
%   outcome         n-length logical vector: 1 if event happened, 0 if not
% Required parameters:
%   prevalence      Prevalence table by age bracket and, optionally, sex
%   default         Prevalence if out of age bracket range
%   ageBracTable    Table of age brackets. 2 columns for lower & upper lims
%   age             n-length vector of ages
%   sex             n-length vector of sexes (1 for male, 2 for female), or
%                   empty vector [] if not sex-dependent
%
% If sex is an empty vector, prevalence must only have 1 column
% If sex is provided, prevalence must have 2 columns (one for each sex)

% If age is empty, then just return an empty set
if isempty(age)
    outcome = false(0,1);
    return
end

% Make sure they provided the right size prevalence table
numAgeBins = size(prevalence,1);
if isempty(sex)
    if size(prevalence,2) ~= 1
        error('simulateEvent: Incorrect prevalence table size (expecting m x 1)');
    end
    sex = 1;
else
    if size(prevalence,2) ~= 2
        error('simulateEvent: Incorrect prevalence table size (expecting m x 2)');
    end
end

% Get number of people
numPpl = length(age);

% Initialize everybody to false
outcome = false(numPpl,1);

% Find age bracket of everybody, and whether they were in range
% ageBrac = ageBracketMaker(ageBracTable,age);
% inRange = ageBrac >= 0;
% numInRange = sum(inRange);
% 
% % Generate linear index
% % (see http://www.mathworks.com/help/techdoc/math/f1-85462.html)
% % to use for looking up prevalence table
% prevalenceBinIndex = ageBrac;
% if (size(prevalence,2) == 2)
%     prevalenceBinIndex(sex==2) = prevalenceBinIndex(sex==2) + numAgeBins;
% end
% 
% % Simulate random event for the people in range
% outcome(inRange) = rand(numInRange,1) < prevalence(prevalenceBinIndex(inRange));
% % Simulate random event for the people out of range
% outcome(~inRange) = rand(numPpl-numInRange,1) < default;

%initialize lookup table
prevLookup = default*ones(110,size(prevalence,2));

%make lookup table
ageBracTable(ageBracTable == Inf) = 109;
for ii = 1:size(ageBracTable,1)
    prevLookup(ageBracTable(ii,1)+1:ageBracTable(ii,2)+1,1) = prevalence(ii,1);
    if size(prevalence,2) > 1
        prevLookup(ageBracTable(ii,1)+1:ageBracTable(ii,2)+1,2) = prevalence(ii,2);
    end
end

% Find outcomes
fate = rand(numPpl,1);
if size(prevalence,2) == 1
    for ii = 1:numPpl
        outcome(ii) = fate(ii) < prevLookup(age(ii)+1,1);
    end
else
    for ii = 1:numPpl
        outcome(ii) = fate(ii) < prevLookup(age(ii)+1,sex(ii));
    end
end
% Done!
