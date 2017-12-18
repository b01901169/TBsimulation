function latToActAnalysis(folderName)
%latToActAnalysis(foldername)
%foldername is a string
%give this function where your 230-year, recordUnit = 1 month file lives 
%and it will throw off some latent to activation stats in a text file in
%that same folder

diary([folderName '/latToActAnalysis.txt'])
diary on

%hard coded for a 230 year run 

%230 year run with no treatment
mainDirectory = pwd;
cd(folderName);

load('recNo.mat');
%recNo = double(intRecNo(:,(75*12):end) );  %truncated for size issues TEMPORARY HARDCODING
recNo = double(intRecNo );  %since recNo is now only written after the burn in period


healthMat = int8(floor( (  mod(recNo ,48) )/6 ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('deleting dead and unborn ppl')
%%%delete all the dead people who are dead or unborn for the analysis period
size(recNo)
recNo2 = zeros(size(recNo,1),size(recNo,2));
nonempty = false(size(recNo,1),1);
for i = 1:length(recNo)
    if (find(healthMat(i,:) == 6,1) ~= 1 ) %not dead to start
        if(  find(healthMat(i,:) == 5, 1,'last') < (50*12) ) %born within the first 50 years
            recNo2(i,:) = recNo(i,:);
            nonempty(i,1) = true;
        end        
    end
end
clear recNo healthMat;
recNo = recNo2(nonempty,:);
clear recNo2 nonemtpy;
disp('size of recNo if keep only 75+ years of nondead and nonunborn records')
size(recNo)

healthMat = int8(floor( (  mod(recNo ,48) )/6 ));
%trtmtVal = int8(mod(recNo,6));
%sexMat = int8((1/48)*( (mod(recNo,144)) - 6*double(healthMat) - double(trtmtVal))); %sex
smokingMat = int8(floor(mod(recNo,43632)/14544));  %smoking
ageMat =  int8( floor( mod(recNo,14544)/144 ) );   %age

%check range of health Mat
disp('age range')
unique(ageMat)

disp('health range')
unique(healthMat)

% disp('treatment range')
% unique(trtmtVal)
% 
% disp('sex range')
% unique(sexMat)

disp('smoking range')
unique(smokingMat)

%make indicator vectors for who was born, had TB, etc
indicatorMat = zeros(size(recNo,1),5);
indicatorMat(:,1) = ((sum((healthMat == 5), 2)) > 0);  %1 if row includes before person was bornm 0 otherwise 
indicatorMat(:,2) = ((sum((healthMat == 6), 2)) > 0);  % if died
indicatorMat(:,3) = ((sum((healthMat == 1 | healthMat == 2), 2)) > 0);  %if had latTB
indicatorMat(:,4) = ((sum((healthMat == 3 | healthMat == 4), 2)) > 0);  % if had actTB
indicatorMat(:,5) = (indicatorMat(:,1) == 1 & indicatorMat(:,2) == 1);  % if has full life history

timingMat = zeros(size(recNo,1),6,'int32');
for i = 1:size(recNo,1)
    %when born
    a = find(healthMat(i,:)==5,1,'last');
    if size(a,2) > 0
        timingMat(i,1) = a; %birth column
    end
    
    %when died
    b = find(healthMat(i,:)==6,1);
    if size(b,2) > 0
        timingMat(i,2) = b; %death column
    end

    %get duration of TB free life for latTB ppl
    c = find(healthMat(i,:)==1 | healthMat(i,:)==2 ,1);
    if size(c,2) > 0
        timingMat(i,3) = c ; %first latTB column
    end
    
    %get duration of latTB life for actTB ppl
    d = find(healthMat(i,:)==3 | healthMat(i,:)==4 ,1);
    if size(d,2) > 0
        timingMat(i,4) = d; %first actTB column
    end
    
%     %get duration of latTB life for actTB ppl
%     e = find(ageMat(i,:)==18,1);
%     if size(e,2) > 0
%         timingMat(i,5) = e; %18 yrs old
%     end
%     
%     f =  find(smokingMat(i,:)==1,1);
%     if size(f,2) > 0
%         timingMat(i,6) = f; %first time smoking
%     end    
end

duration = zeros(size(recNo,1),5);
duration(:,1) = timingMat(:,2) - timingMat(:,1);  %total lifespan in months
duration(:,2) = timingMat(:,3) - timingMat(:,1);  %time healthy
duration(:,3) = timingMat(:,4) - timingMat(:,3);  %time latent
duration(:,4) = timingMat(:,2) - timingMat(:,4);  %time active
duration(:,5) = timingMat(:,2) - timingMat(:,5);  %duration of life given lived to age 18

%make cohort
cohort = (timingMat(:,1) >= 1 & timingMat(:,1) <= (5*12) & timingMat(:,2) > 0 );  %HARD CODING for a 230 year run (take ppl born in 1st yr to 5th year)
smokers = ((sum(smokingMat, 2)) > 0);

%============NonTB life expectancy=================

noTBcohort = (timingMat(:,1) >= 1 & timingMat(:,1) <= (5*12) & timingMat(:,2) > 0 & timingMat(:,4) == 0);  %HARD CODING

ActDuration = duration(noTBcohort,1);
avelife = mean(ActDuration);
medlife = median(ActDuration);
disp('average lifespan in recordUnits:')
avelife
disp('median lifespan in recordUnits:')
medlife
disp('for N over:')
sum(noTBcohort)

%do male and female

% %WTF WHY DOES THIS NOT WORK WTF WTF
% %do smoking and nonsmoking
% %need to make a better cohort: 18-yr-old smokers and 18-yr-old nonsmokers
% smokingB4 = false(size(recNo,1), 1);
% smokedB4 = (timingMat(:,6) < timingMat(:,5));
% 
% noTBSmokingCohort = (timingMat(:,1) > (1) & timingMat(:,1) <= (50*12) & timingMat(:,2) > 0 & duration(:,5) > 0 & timingMat(:,4) == 0 & cohort(:,1) == 1 & smokedB4 ==1);  %HARD CODING
% ActDuration = duration(noTBSmokingCohort,5);
% avelife = mean(ActDuration);
% disp('average lifespan in recordUnits:')
% avelife
% 
% %live until at least 18 and never smoke
% noTBNonsmokingCohort = (timingMat(:,1) > (1) & timingMat(:,1) <= (50*12) & timingMat(:,2) > 0 & duration(:,5) > 0 & timingMat(:,4) == 0 & smokers(:,1) == 0 & cohort(:,1) == 1);  %HARD CODING
% ActDuration = duration(noTBNonsmokingCohort,5);
% avelife = mean(ActDuration);
% disp('average lifespan in recordUnits:')
% avelife

%============TB ppl=================

latTBdeathKeepers = (indicatorMat(:,3) == 1 & cohort(:,1) == 1);
actTBdeathKeepers = (indicatorMat(:,4) == 1 & cohort(:,1) == 1 );

%percentage activate
numHadLatent = sum(indicatorMat(latTBdeathKeepers,3));  %everyone who had latent and died
numHadActive = sum(indicatorMat(latTBdeathKeepers,4)); 
percActivate = numHadActive/numHadLatent ;
disp('all TB ppl: percentage activated:')
percActivate
disp('for latent pop of:')
sum(latTBdeathKeepers)

%duration of life post actTB
actTBpplActDuration = duration(actTBdeathKeepers,4);
aveActTBlife = mean(actTBpplActDuration);
disp('all TB ppl: average lifespan post activation in recordUnits:')
aveActTBlife
disp('for act pop of:')
sum(actTBdeathKeepers)


%=====stratification by smoking
latTBdeathSmokers = (indicatorMat(:,3) == 1 & smokers(:,1) == 1 & cohort(:,1) == 1);  %check vector lengths
actTBdeathSmokers = (indicatorMat(:,4) == 1 & smokers(:,1) == 1 & cohort(:,1) == 1);

numHadLatentSmokers = sum(indicatorMat(latTBdeathSmokers,3));  %everyone who had latent and died
numHadActiveSmokers = sum(indicatorMat(latTBdeathSmokers,4)); 
percActivateSmokers = numHadActiveSmokers/numHadLatentSmokers ;
disp('smoking TB ppl: percentage activated smokers:')
percActivateSmokers
disp('for latent smoker pop of:')
sum(latTBdeathSmokers)


actTBpplActDuration = duration(actTBdeathSmokers,4);  %duration of life post actTB
aveActTBlifeSmokers = mean(actTBpplActDuration);
disp('smoking TB ppl: average lifespan smokers post activation in recordUnits:')
aveActTBlifeSmokers
disp('for act smoker pop of:')
sum(actTBdeathSmokers)

latTBdeathNonsmokers = (indicatorMat(:,3) == 1 & smokers(:,1) == 0 & cohort(:,1) == 1);
actTBdeathNonsmokers = (indicatorMat(:,4) == 1 & smokers(:,1) == 0 & cohort(:,1) == 1);

numHadLatentNonsmokers = sum(indicatorMat(latTBdeathNonsmokers,3));  %everyone who had latent and died
numHadActiveNonsmokers = sum(indicatorMat(latTBdeathNonsmokers,4)); 
percActivateNonsmokers = numHadActiveNonsmokers/numHadLatentNonsmokers;
disp('nonsmoking TB ppl: percentage activated smokers:')
percActivateNonsmokers
disp('for latent nonsmoker pop of:')
sum(latTBdeathNonsmokers)

actTBpplActDuration = duration(actTBdeathNonsmokers,4);  %duration of life post actTB
aveActTBlifeNonsmokers = mean(actTBpplActDuration);
disp('nonsmoking TB ppl: average lifespan smokers post activation in recordUnits:')
aveActTBlifeNonsmokers
disp('for act nonsmoker pop of:')
sum(actTBdeathNonsmokers)


% 
% %finding force of infection
% %count how many people in each age grp at last time period
% ageIndicator = false(size(ageMat,1),15);
% 
%     alive = (healthMat(:,end) ~=5 & healthMat(:,end) ~=6);
%     lastAge = ageMat(:,end);
%     ageIndicator(:,1) = (lastAge <= 4 & alive & cohort);
%     ageIndicator(:,2) = (lastAge >= 5 & lastAge <= 9 & alive & cohort);
%     ageIndicator(:,3) = (lastAge >= 10 & lastAge <= 14 & alive & cohort);
%     ageIndicator(:,4) = (lastAge >= 15 & lastAge <= 19 & alive & cohort);
%     ageIndicator(:,5) = (lastAge >= 20 & lastAge <= 24 & alive & cohort);
%     ageIndicator(:,6) = (lastAge >= 25 & lastAge <= 29 & alive & cohort);
%     ageIndicator(:,7) = (lastAge >= 30 & lastAge <= 34 & alive & cohort);
%     ageIndicator(:,8) = (lastAge >= 35 & lastAge <= 39 & alive & cohort);
%     ageIndicator(:,9) = (lastAge >= 40 & lastAge <= 44 & alive & cohort);
%     ageIndicator(:,10) = (lastAge >= 45 & lastAge <= 49 & alive & cohort);
%     ageIndicator(:,11) = (lastAge >= 50 & lastAge <= 54 & alive & cohort);
%     ageIndicator(:,12) = (lastAge >= 55 & lastAge <= 59 & alive & cohort);
%     ageIndicator(:,13) = (lastAge >= 60 & lastAge <= 64 & alive & cohort);
%     ageIndicator(:,14) = (lastAge >= 65 & lastAge <= 69 & alive & cohort);
%     ageIndicator(:,15) = (lastAge >= 70 & alive & cohort);
%     
%     numPplInAge = sum(ageIndicator);
%     fracPplInAge = numPplInAge ./[sum(alive), sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive),sum(alive) ]
%     
%     disp('histogram: frequency of people in each age, with bins = 0:5:100')
%     histc(ageMat(:,120),[0:5:100])'
%     
%     contactsPerDay = (fracPplInAge)*[...
%         5.45
%         7.64
%         6.66
%         6.44
%         5.71
%         5.74
%         6.02
%         6.35
%         5.66
%         4.42
%         3.31
%         3.29
%         2.61
%         1.95
%         2.84
%         ];
% 
% forceOfInfection = contactsPerDay*365*0.0029;
% 
% disp('The force of infection is (for transmission scaling of 0.0029)')
% forceOfInfection

%%%%%plotter%%%
% freqCounter = zeros(max(ActDuration),1);
% for i = 1:size(ActDuration,1)
%     freqCounter(ActDuration(i,1),1) = freqCounter(ActDuration(i,1),1)  +1;
% end
% plot(freqCounter)

diary off

cd(mainDirectory);

