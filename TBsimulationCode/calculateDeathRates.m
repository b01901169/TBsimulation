function [chances, monthlyDeathProb, yearlyDeathRate, fiveYearDeathProb, yearlyDeathRateSmokingAvg] = ...
    calculateDeathRates(aliveMat, deathMat)

% [chances, 

%  monthlyDeathProb, 

%  yearlyDeathRate, 

%  fiveYearDeathProb,

%  yearlyDeathRateSmokingAvg 

% ] = calculateDeathRates(aliveMat, deathMat)

%

% Calculates death rates and just displays them (modify this

% function if you want to do anything with them)



% Sum over time to get the total number of times that a person of a certain

% age has died due to natural causes, as well as the number of

% opportunities (i.e. months) that people of a given age have had to die

deaths  = sum(deathMat,1);

chances = sum(aliveMat,1);

% Now get rid of the first dimension

deaths  = squeeze(deaths);

chances = squeeze(chances);

% And divide to get monthly probabilities.

monthlyDeathProb = deaths ./ chances;

% Convert to yearly death rate

yearlyDeathRate = -log(1-monthlyDeathProb) * 12;

% Also convert to 5-year death prob

fiveYearDeathProb = 1 - (1-monthlyDeathProb).^12;

% Display in a table

smokingTitles = {'Nonsmokers','Smokers'};

for smoke = 1:2

    for sex = 1:2

        fprintf('%s\n',smokingTitles{smoke});

        fprintf('Sex = %d\n',sex);

        fprintf('Age | Pop*months | Month prob | Year rate | 5-year prob \n');

        fprintf('----+------------+------------+-----------+-------------\n');

        for ii = 1:101

            fprintf('%3d | %10d | %010.8f | %09.7f | %09.7f\n',...
                ii-1, chances(ii,sex,smoke), monthlyDeathProb(ii,sex,smoke), ...
                yearlyDeathRate(ii,sex,smoke), fiveYearDeathProb(ii,sex,smoke));

        end

        fprintf('\n');

    end

    fprintf('--------------------------------------------------------\n');

end



% Do the smoking average

monthlyDeathProbSmokingAvg = sum(deaths,3) ./ sum(chances,3);

yearlyDeathRateSmokingAvg = -log(1-monthlyDeathProbSmokingAvg) * 12;



