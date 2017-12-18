function plotByAge(stateMat, Page, toPlot)
% plotPrevByAge(ppl, Page, toPlot)
%
% Plots a plot of people
%   ppl     people matrix
%   Page    Index of age
%   toPlot  Indices to plot in people matrix

% Things to plot
numStats = length(toPlot);

% Correction factor
correctionFactor = min(stateMat(:,toPlot));

% Age bins
ageLowerLim = (0:5:95)';
ageUpperLim = (5:5:100)';
ageBracTable = [ageLowerLim, ageUpperLim];
numBins = length(ageBracTable);

% Generate stats
numPpl = length(stateMat);
stats = zeros(numBins,numStats + 1);
ageBrac = ageBracketMaker(ageBracTable,stateMat(:,Page));
for ageIndex = 1:numBins
    % Number of people in that age bracket that have a specified stat
    for statIndex = 1:numStats
        stats(ageIndex,statIndex) = sum(stateMat(ageBrac == ageIndex,toPlot(statIndex)) - correctionFactor(statIndex));
    end
    % Total number of people in that age bracket
    stats(ageIndex,numStats+1) = sum(ageBrac == ageIndex);
end

warning off   %turn divide by zero warning off

% Plot
figure
subplot(2,1,1)
plot(mean(ageBracTable,2), stats);
leg = legend(num2str(reshape(toPlot,numStats,1)));
set(leg,'Location','NorthEastOutside');
xlabel('age');
ylabel('no. ppl');
subplot(2,1,2)
plot(mean(ageBracTable,2), stats(:,1:numStats) ./ (stats(:,end) * ones(1,numStats)));
xlabel('age');
ylabel('% ppl in age group')

warning on   %turn divide by zero warning on