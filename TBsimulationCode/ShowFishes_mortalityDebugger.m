function ShowFishes_mortalityDebugger(fishFolder, runYear)

% Script to plot supposed death rates vs. that made by simulation



folderName = fishFolder;

plotFixerFold = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs'; %this is where plotfixer lives;



ShowFishes_whoData;



% Calculate from simulation output

currentDirectory = pwd;

cd(fishFolder);

load fishMatrices.mat

cd(currentDirectory);

[chances, monthlyDeathProb, yearlyDeathRate, fiveYearDeathProb, yearlyDeathRateSmokingAvg] = ...
    calculateDeathRates(aliveAgesMat, deathAgesMat);



%% Plot yearly death rates

figure(1)

% Males

subplot(2,1,1);

sex = 1;

whoX = reshape(whoAgeBrac',2*length(whoAgeBrac),1);

whoY = reshape(repmat(who1990YearlyDeathRate(:,sex),[1 2])',2*length(whoAgeBrac),1);

plot(whoX,whoY,'k-');

hold on

smoke = 1; % smoke = 1 is nonsmoking, 2 is smoking

plot(0:100,yearlyDeathRate(:,sex,smoke),'c.')

hold on

smoke = 2; % smoke = 1 is nonsmoking, 2 is smoking

plot(0:100,yearlyDeathRate(:,sex,smoke),'r.')

hold on

plot(0:100,yearlyDeathRateSmokingAvg(:,sex), 'b.')

hold off

legend('WHO','Simulation Nonsmokers','Simulation Smokers','Simulation AveSmoking','Location','North');

xlabel('Age');

ylabel('Yearly death rate');

title('Males, using 1990 data');

% Females

subplot(2,1,2);

sex = 2;

whoX = reshape(whoAgeBrac',2*length(whoAgeBrac),1);

whoY = reshape(repmat(who1990YearlyDeathRate(:,sex),[1 2])',2*length(whoAgeBrac),1);

plot(whoX,whoY,'k-');

hold on

smoke = 1;

plot(0:100,yearlyDeathRate(:,sex,smoke),'c.')

hold on

smoke = 2;

plot(0:100,yearlyDeathRate(:,sex,smoke),'r.')

hold on

plot(0:100,yearlyDeathRateSmokingAvg(:,sex), 'b.')

hold off

xlabel('Age');

ylabel('Yearly death rate');

title('Females, using 1990 data');

cd(plotFixerFold);

plotfixer;

cd(currentDirectory);

print('-dpng','-r70',[folderName '/yearlyDeathRate.png']);



%% Wow, that was really close, let's look at it on a log plot

figure(2)

% Males

subplot(2,1,1);

sex = 1;

whoX = reshape(whoAgeBrac',2*length(whoAgeBrac),1);

whoY = reshape(repmat(who1990YearlyDeathRate(:,sex),[1 2])',2*length(whoAgeBrac),1);

semilogy(whoX,whoY,'k-');

hold on

smoke = 1; % smoke = 1 is nonsmoking, 2 is smoking

semilogy(0:100,yearlyDeathRate(:,sex,smoke),'c.')

hold on

smoke = 2; % smoke = 1 is nonsmoking, 2 is smoking

semilogy(0:100,yearlyDeathRate(:,sex,smoke),'r.')

hold on

semilogy(0:100,yearlyDeathRateSmokingAvg(:,sex), 'b.')

hold off

legend('WHO','Simulation Nonsmokers','Simulation Smokers','Simulation AveSmoking','Location','North');

xlabel('Age');

ylabel('Yearly death rate');

title('Males, using 1990 data');

% Females

subplot(2,1,2);

sex = 2;

whoX = reshape(whoAgeBrac',2*length(whoAgeBrac),1);

whoY = reshape(repmat(who1990YearlyDeathRate(:,sex),[1 2])',2*length(whoAgeBrac),1);

semilogy(whoX,whoY,'k-');

hold on

smoke = 1;

semilogy(0:100,yearlyDeathRate(:,sex,smoke),'c.')

hold on

smoke = 2;

semilogy(0:100,yearlyDeathRate(:,sex,smoke),'r.')

hold on

semilogy(0:100,yearlyDeathRateSmokingAvg(:,sex), 'b.')

hold off

xlabel('Age');

ylabel('Yearly death rate');

title('Females, using 1990 data');

cd(plotFixerFold);

plotfixer;

cd(currentDirectory);

print('-dpng','-r70',[folderName '/yearlyDeathRate_logPlot.png']);



%% Extra fishy plot

figure(3)

load fishMatrices_extraFishy.mat

imagesc(A);

axis image





%% Figure out ages at death

deathsByAge = squeeze(sum(deathAgesMat,1));

% Norm it to per 100,000

deathsNormed = zeros(size(deathsByAge));

deathsNormed_AveSmoking = zeros(size(deathsByAge,1), size(deathsByAge,2));

for sex = 1:2

    for smoke = 1:2

        deathsNormed(:,sex,smoke) = ...
            100000 * deathsByAge(:,sex,smoke) / sum(deathsByAge(:,sex,smoke));

    end

    deathsNormed_AveSmoking(:,sex) = 100000 * sum(deathsByAge(:,sex,:),3) / sum(  sum(deathsByAge(:,sex,:),3) );

end

% Convert into WHO age bracket bins

deathsByWHOAgeBrac = zeros(length(whoAgeBrac),2,2);

deathsByWHOAgeBrac_AveSmoking = zeros(length(whoAgeBrac),2);

for age = 0:100

    ageBrac = find(whoAgeBrac(:,1) <= age & whoAgeBrac(:,2) >= age);

    deathsByWHOAgeBrac(ageBrac,:,:) = ...
        deathsByWHOAgeBrac(ageBrac,:,:) + deathsNormed(age+1,:,:);

    deathsByWHOAgeBrac_AveSmoking(ageBrac,:) = deathsByWHOAgeBrac_AveSmoking(ageBrac,:) + deathsNormed_AveSmoking(age+1,:);

end



% Use the centers of the age brackets for the x-axis when plotting

ageBracCenters = mean(whoAgeBrac,2);



% Plot

figure(4)

% males

subplot(2,1,1)

sex = 1;

plot(ageBracCenters, who1990deathsByAge(:,sex), 'ko', ...
    ageBracCenters, deathsByWHOAgeBrac(:,sex,1), 'c.',...
    ageBracCenters, deathsByWHOAgeBrac(:,sex,2), 'r.',...
    ageBracCenters, deathsByWHOAgeBrac_AveSmoking(:,sex), 'b.' );  %third index on deathsByWHOAgeBrac is smoking.  1 is nonsmoke, 2 is smoke.

xlabel('Age');

ylabel('#/100k that die in that age bracket');

legend('WHO','Simulation Nonsmokers','Simulation Smokers','Simulation AveSmoking','Location','NorthWest');

title('Males');

% females

subplot(2,1,2)

sex = 2;

plot(ageBracCenters, who1990deathsByAge(:,sex), 'ko', ...
    ageBracCenters, deathsByWHOAgeBrac(:,sex,1), 'c.',...
    ageBracCenters, deathsByWHOAgeBrac(:,sex,2), 'r.',...
    ageBracCenters, deathsByWHOAgeBrac_AveSmoking(:,sex), 'b.');  %third index on deathsByWHOAgeBrac is smoking.  1 is nonsmoke, 2 is smoke.

xlabel('Age');

ylabel('#/100k that die in that age bracket');

title('Females');

% Annotations

% for sex = 1:2

%     subplot(2,1,sex);

%     text(78,1.8e4,'Life Exp. from 15');

%     x = 93; y = 1.5e4; dy = -.25e4;

%     text(x,y,sprintf('%.1f',LE_WHO(sex)),'Color','k'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_5yr(sex,1)),'Color','c'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_5yr(sex,2)),'Color','r'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_5yr_AveSmoking(sex)),'Color','b'); y = y + dy;

% end

cd(plotFixerFold);

plotfixer;

cd(currentDirectory);

print('-dpng','-r70',[folderName '/numDeadbyAge.png']); 



%% Also do a plot by year (instead of by age bracket) because that is also

% interesting to look at



% Plot

figure(5)

for sex = 1:2

    subplot(2,1,sex);

    whoX = reshape(whoAgeBrac',2*length(whoAgeBrac),1);

    whoY = reshape(repmat(who1990deathsByAge(:,sex),[1 2])',2*length(whoAgeBrac),1);

    % Need to scale this by the number of years in each age bracket

    % Infant mortality gets scaled by 4 to keep it from stretching out the y-scale

    % I feel kinda bad about hardcoding this...

    whoY(1:4) = whoY(1:4) / 4;

    whoY(5:end) = whoY(5:end) / 5;

    % Plot

    plot(whoX,whoY,'k-');

    hold on

    plot(1:100,deathsNormed(2:end,sex,1),'c.')  %1 is nonSmoking

    plot(1:100,deathsNormed(2:end,sex,2),'r.')  %2 is Smoking

    plot(1:100,deathsNormed_AveSmoking(2:end,sex),'b.')  %2 is Smoking

    plot(0,deathsNormed(1,sex,1)/4,'c.')        %1 is nonSmoking

    plot(0,deathsNormed(1,sex,2)/4,'r.')        %2 is Smoking

    plot(0,deathsNormed_AveSmoking(1,sex)/4,'b.')        %2 is Smoking

    hold off

    xlabel('Age');

    ylabel('#/100k that die at that age');

    title(sprintf('Sex = %d',sex));

    % Life expectancy

%     text(80,3500,'Life Exp. from 15');

%     x = 93; y = 3000; dy = -500;

%     text(x,y,sprintf('%.1f',LE_WHO(sex)),'Color','k'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_1yr(sex,1)),'Color','c'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_1yr(sex,2)),'Color','r'); y = y + dy;

%     text(x,y,sprintf('%.1f',LE_1yr_AveSmoking(sex)),'Color','b'); y = y + dy;

end

subplot(2,1,1);

legend('WHO','Simulation Nonsmokers','Simulation Smokers','Simulation AveSmoking','Location','NorthWest');

subplot(2,1,2);

text(4,7000,'Divided infant mortality by 4 to keep scale manageable');

cd(plotFixerFold);

plotfixer;

cd(currentDirectory);

print('-dpng','-r70',[folderName '/numDeadbyYear.png']);



close all







%%%%%%%%%%%%%Calculating Life Expectancy %%%%%%%%%%%%%%%%%%%%%%%%

LE_WHO = zeros(2,1);

LE_WHO_startAt30 = zeros(2,1);



%make age brackets for life expectancy

ages = ageBracCenters(1:end);  % baseline LE calculation age brackets

ages_from30 = ageBracCenters(8:end);  % Don't include the infant mortality - start at 30 (age bracket 8)



%find WHO life expectancy

for sex = 1:2

    LE_WHO(sex) = calculateLifeExpectancy(ages,...
        who1990deathsByAge(1:end,sex));

    LE_WHO_startAt30(sex) = calculateLifeExpectancy(ages_from30,...
        who1990deathsByAge(8:end,sex));

end



% First, calculate life expectancy

% Start at 15

ageCenters = (0:100)' + 0.5;

ageCenters_from30 = (30:100)' + 0.5;

LE_1yr = zeros(2,2);

LE_1yr_AveSmoking = zeros(2,1);

for sex = 1:2

    for smoke = 1:2

        LE_1yr(sex,smoke) = calculateLifeExpectancy(ageCenters_from30, ...
            deathsNormed(30+1:end,sex,smoke));  %for smoking, start LE calculation at age 30

    end

    LE_1yr_AveSmoking(sex) = calculateLifeExpectancy(ageCenters, ...
        deathsNormed_AveSmoking(0+1:end,sex)/100000);

end



LE_output = [LE_WHO(1),LE_1yr_AveSmoking(1),LE_WHO_startAt30(1),LE_1yr(1,1) ,LE_1yr(1,2) ; LE_WHO(2),LE_1yr_AveSmoking(2),LE_WHO_startAt30(2),LE_1yr(2,1) ,LE_1yr(2,2)];

tablePrinter('WHO LE, Simu LE, WHO LE from 30, nonsmokers Simu LE from 30, smokers Simu LE from 30, row1 male and row2 female', LE_output, 'LifeExp', folderName);





%%%%%%%%%%%%%Calculating Life Expectancy After TB activates%%%%%%%%%%%%%%%%%%%%%%%%



%gotta chantge this



%find life expectancy

alive = squeeze(sum(aliveActTBAgesMat,1));

alive_after = squeeze(sum(aliveActTBAgesMat_after,1));

died = zeros(size(alive));



%got rid of the timePeriod dimension.  now have: dim 1 = mons since activate,

%dim 2 = nonsmoker/smoker, dim3 = under 30/over30



for col = 1:size(alive,1)-1  %minus 1 since aliveActTBAgesMat is truncated and ppl with 250 months since activation stay there until they die (counted repeatedly)

    died(col,:) = alive(col,:) - alive_after(col,:);  

end

deathFrac = sum(sum(died,2),3) ./ sum(sum(sum(died,2),1),3);



diedOver30 = died(:,:,2);

nonSmokeDeathFrac = diedOver30(:,1,:) ./ sum(diedOver30(:,1,:),1);

smokeDeathFrac = diedOver30(:,2,:) ./ sum(diedOver30(:,2,:),1);



deathFrac(isnan(deathFrac)) = 0 ;

nonSmokeDeathFrac(isnan(nonSmokeDeathFrac)) = 0 ;

smokeDeathFrac(isnan(smokeDeathFrac)) = 0 ;



generalLifeEx = [1:size(died,1)] * ( deathFrac );

nonsmokersLifeEx = [1:size(died,1)] * (  nonSmokeDeathFrac  );

smokersLifeEx  = [1:size(died,1)] * (  smokeDeathFrac );



%find the percentage of people who activate

percActivate = size(find(MonSinceAct~=0))/size(find(MonSinceLat~=0)) ;



tablePrinter('general, nonsmokers, smokers, percActivate, totNumPpl', [generalLifeEx,nonsmokersLifeEx,smokersLifeEx,percActivate,sum(sum(sum(died,2),1),3)], 'postActivationLifeExp', folderName);



nonZeroMonSinceActTruncations = size(find(MonSinceActTruncations~=0));



tablePrinter('num of nonZero elements MonSinceActTruncations', nonZeroMonSinceActTruncations, 'nonZeroMonSinceActTruncations', folderName);

