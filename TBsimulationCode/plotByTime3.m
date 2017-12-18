function output  = plotByTime3(recNo, subsetFlag, subsetIndicator, unit, ageFlag)

%function plotByTime3(recNo, subsetFlag, subsetIndicator, unit)

%Plots the frequency of ACTIVE TB outcomes at every time period.

%

%RECNO is the compressed stateMatrix with all the information in it.

%SUBSETFLAG is 0 or 1.  0 is not need to separate by sex, smoking, urb/rur.

%subsetIndicator is ignored if subset == 0.  Otherwise, flag = 1 for sex, 2 for

%smoking, 3 for urbanRural.  

%

%UNIT is number of months between bars on graph (ex, unit = 12 is time is in years)

%ageFlag is 0 if the y-axis should be treatment/health, 1 if the y-axis

%should be age

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%tic

%fprintf('  Plotting (%d,%d,%d,%d)\n',subsetFlag,subsetIndicator,unit,ageFlag);



%initialize vars

numPpl = size(recNo,1);

totPeriods = size(recNo,2);

numCodes = max(max(recNo));

graphMat = [];

propMat = [];



%make the definition of the categories

codeNumMat = zeros(numCodes,1);

%make matrix of code values

for i = 1: numCodes

    codeNumMat(i,1) = i;

end



%always need a healthMat

healthMat = int8(floor( (  mod(codeNumMat,48) )/6 ));



%make y-axis variable

if (ageFlag == 0)

    %make indicies easier to read;

    numHealthOutcomes = 13;

    healthy = 1; latentSens = 2; latentMdr = 3; activeSens =4; activeMdr = 5; notYetBorn = 6; dead = 7;

    catI_sens = 8; catII_sens = 9; catIII = 10; catIV = 11; catI_MDR = 12; catII_MDR = 13;



    toPlot = [activeSens; catI_sens; catII_sens; activeMdr; catI_MDR; catII_MDR];

    legendText = {'actSens', 'catI sens', 'catII sens', 'actMdr', 'catI MDR',  'catII MDR'};

    header1 = ['actSens , catI sens , catII sens , actMdr , catI MDR , catII MDR' ]; 

    header2 = [',']; %no proportions for plotByTime3

    yVarSt = 'Frequencies of ActTB Health Outcomes';

    

    %decode the codes

    trtmtVal = int8(mod(codeNumMat,6));



    % Generate healthTreatVal

    healthTreatVal = zeros(numCodes, 1, numHealthOutcomes);

    healthTreatVal(:, 1,activeSens) = (healthMat ==3 & trtmtVal == 0 );

    healthTreatVal(:, 1,activeMdr) = (healthMat == 4 & trtmtVal == 0 );

    

    healthTreatVal(:, 1,catI_sens) = (trtmtVal == 1 & healthMat == 3 );

    healthTreatVal(:, 1,catII_sens) = (trtmtVal == 2 & healthMat == 3  );

    healthTreatVal(:, 1,catI_MDR) = (trtmtVal == 1 & healthMat == 4  );

    healthTreatVal(:, 1,catII_MDR) = (trtmtVal == 2 & healthMat == 4  );



    healthTreatVal(:, 1,healthy) = (healthMat == 0);

    healthTreatVal(:, 1,latentSens) = (healthMat == 1 );

    healthTreatVal(:, 1,latentMdr) = (healthMat == 2);

    healthTreatVal(:, 1,catIII) = (trtmtVal == 3 ); %all catIII

    healthTreatVal(:, 1,catIV) = (trtmtVal == 4  );  %all catIV

    healthTreatVal(:, 1,notYetBorn) = (healthMat == 5);

    healthTreatVal(:, 1,dead) = (healthMat == 6);

    

    

    

else  %ageFlag == 1, and want y-axis to be age

    %make indicies easier to read;

    numHealthOutcomes = 15;

    sens_under20 = 1; sens_From20to40 = 2; sens_From40to60 =3; sens_From60to80 =4; sens_From80to100=5;...
        MDR_under20=6; MDR_From20to40 =7 ; MDR_From40to60 = 8; MDR_From60to80 = 9; MDR_From80to100 = 10;...
        notYetBorn = 11; dead = 12; healthy = 13; latentSens = 14; latentMdr = 15;



    toPlot = [sens_under20; sens_From20to40; sens_From40to60; sens_From60to80; sens_From80to100;...
        MDR_under20; MDR_From20to40; MDR_From40to60; MDR_From60to80; MDR_From80to100];

    legendText = {'sens under20',' sens From20to40',' sens From40to60',' sens From60to80',' sens From80to100',...
        'MDR under20',' MDR From20to40',' MDR From40to60',' MDR From60to80',' MDR From80to100'};

    header1 = ['sens_under20',' sens_From20to40',' sens_From40to60',' sens_From60to80',' sens_From80to100',...
        'MDR_under20',' MDR_From20to40',' MDR_From40to60',' MDR_From60to80',' MDR_From80to100' ]; 

    header2 = [',']; %no proportions for plotByTime3



    yVarSt = 'Age Structure';



    ageMat = int8( floor( mod(codeNumMat,14544)/144 ) );



    healthTreatVal = zeros(numCodes, 1, numHealthOutcomes);

    healthTreatVal(:,1,sens_under20) = (ageMat <= 20 & healthMat == 3 );

    healthTreatVal(:,1,sens_From20to40) = (ageMat > 20 & ageMat <= 40 & healthMat == 3 );

    healthTreatVal(:,1,sens_From40to60) = (ageMat > 40 & ageMat <= 60 & healthMat == 3 );

    healthTreatVal(:,1,sens_From60to80) = (ageMat > 60 & ageMat <= 80 & healthMat == 3  );

    healthTreatVal(:,1,sens_From80to100) = (ageMat > 80 & ageMat <= 100 & healthMat == 3  );

    

    healthTreatVal(:,1,MDR_under20) = (ageMat <= 20 & healthMat ~= 6 & healthMat == 4 );

    healthTreatVal(:,1,MDR_From20to40) = (ageMat > 20 & ageMat <= 40 & healthMat == 4 );

    healthTreatVal(:,1,MDR_From40to60) = (ageMat > 40 & ageMat <= 60 & healthMat == 4 );

    healthTreatVal(:,1,MDR_From60to80) = (ageMat > 60 & ageMat <= 80 & healthMat == 4  );

    healthTreatVal(:,1,MDR_From80to100) = (ageMat > 80 & ageMat <= 100 & healthMat == 4  );

    

    healthTreatVal(:, 1,notYetBorn) = (healthMat == 5);

    healthTreatVal(:, 1,dead) = (healthMat == 6);    

    healthTreatVal(:, 1,healthy) = (healthMat == 0);

    healthTreatVal(:, 1,latentSens) = (healthMat == 1);

    healthTreatVal(:, 1,latentMdr) = (healthMat == 2);

    

end

%toc



%make subset matracies if necessary

if (subsetFlag == 1)  %with subsetting

    if subsetIndicator == 1  %sex

        subsetMat = int8((1/48)*( (mod(codeNumMat,144)) - 6*double(healthMat) - double(trtmtVal)));

        titleSt = sprintf('%s by %s (lefthand graphs male, righthand female)' ,yVarSt, 'Sex');

        categories = {'male','female'};

    elseif subsetIndicator == 2  %smoking

        subsetMat = int8(floor(mod(codeNumMat,43632)/14544));

        titleSt = sprintf('%s by %s (lefthand graphs nonsmoking, righthand smoking)',yVarSt,'Smoking Status');

        categories = {'nonsmoking','smoking'};

    elseif subsetIndicator == 3 %urban rural

        subsetMat = int8(floor(codeNumMat/43632));

        titleSt = sprintf('%s by %s (lefthand graphs urban, righthand rural)',yVarSt,'Urban/Rural');

        categories = {'urban','rural'};

    else

        disp('no such subset flag in graphing func')

        titleSt = ' ';

        categories = {'',''};

    end

end





%TIME LOOP

for timePeriod = 1: unit : totPeriods

    %reset vars

    codeFreqMat = zeros(numCodes,1);

    proportions = zeros(numCodes,1);



    % sum recNo codes

    for person = 1 : numPpl

        code = recNo(person,timePeriod);

        codeFreqMat(code) = codeFreqMat(code) + 1;

    end

        

    %make freqencies

    if (subsetFlag == 0)  %with no subsetting



        for i = 1:numHealthOutcomes

            healthMarker = (healthTreatVal(:,1,i) == 1);

            healthFreq(1,i) = sum(codeFreqMat(healthMarker,1));

        end

        titleSt = yVarSt;

        

        proportions = double(healthFreq./(numPpl-healthFreq(notYetBorn)-...
            healthFreq(dead) - healthFreq(healthy) - ...
            healthFreq(latentSens) - healthFreq(latentMdr) ));



    elseif (subsetFlag == 1)  %with subsetting



        for j = 1:2         %loop over subset types

            for i = 1:numHealthOutcomes         %loop over health outcomes

                subsetVal = j;

                if subsetIndicator == 2  %if smoking is the subset desired

                    subsetVal = j - 1;

                end

                healthMarker = (healthTreatVal(:,1,i) == 1 & subsetMat == subsetVal);                

                healthFreq(j,i) = sum(codeFreqMat(healthMarker,1));

            end    %end loop over health outcomes

        end  %end loop over subset types

        

        %add empty row

        healthFreq(3,:) = zeros(1, numHealthOutcomes);

        numAlive =  (sum(healthFreq,2)-healthFreq(:,notYetBorn)-...
            healthFreq(:,dead) - healthFreq(:,healthy) - ...
            healthFreq(:,latentSens) - healthFreq(:,latentMdr) ); %this is actually the number of total activeTB

        

        %make the matrix of [numAlive, numAlive,...] for making proportions

        numAliveRepeated = numAlive;

        for repeats = 1 : numHealthOutcomes - 1

            numAliveRepeated = [numAliveRepeated, numAlive];

        end

        

        proportions = double(healthFreq./ numAliveRepeated);    

            %there should be as many columns as there are in healthFreq

            

    end   %end if loop



    %add this time period to the graphMat

    graphMat = [graphMat; healthFreq];

    propMat = [propMat; proportions];

end  %end timePeriod loop



%add an extra empty entry at the end for space

if (size(graphMat,1) == 1)

    graphMat = [graphMat; zeros(length(healthFreq))];

    propMat = [propMat; zeros(length(proportions))];

end

%toc



graphTable = graphMat(:,toPlot');



if (subsetFlag == 0)  %with no subsetting

    graphTable = graphMat(:,toPlot');

    

    %plot the graphMat

    subplot(2,1,1);

    bar(graphTable,'stack')

    leg = legend(legendText);

    set(leg,'Location','NorthEastOutside');

    title(titleSt);

    hold on;



    %plot the propMat

    subplot(2,1,2);

    bar(propMat(:,toPlot'),'stack');

    ylim([0,1]);

    title('(by proportion of alive)')

    

elseif (subsetFlag == 1)  %with subsetting

    %separate the graph so it is nonsmoking/smoking, male/female, etc

    graphMatSubsetVal1 = zeros(size(graphMat,1)/3,size(graphMat,2));

    graphMatSubsetVal2 = zeros(size(graphMat,1)/3,size(graphMat,2));

    propMatSubsetVal1 = zeros(size(propMat,1)/3,size(propMat,2));

    propMatSubsetVal2 =zeros(size(propMat,1)/3,size(propMat,2));

    for row = 1: size(graphMat,1)

        rowDouble = double(row);

        if mod(row,3) == 1

            graphMatSubsetVal1(floor(rowDouble/3)+1,:) = graphMat(row,:);        

            propMatSubsetVal1(floor(rowDouble/3)+1,:) = propMat(row,:);

        elseif mod(row,3) == 2

            propMatSubsetVal2(floor(rowDouble/3)+1,:) = graphMat(row,:);

            graphMatSubsetVal2(floor(rowDouble/3)+1,:) = propMat(row,:);

        end

    end    

    

    %plot the graphMatSubsetVal1

    graphsForPlotting = {graphMatSubsetVal1, graphMatSubsetVal2};

    propsForPlotting = {propMatSubsetVal1, propMatSubsetVal2};

    

    for subTypes = 1:2

        graphTable = graphsForPlotting{subTypes};

        propMat = propsForPlotting{subTypes};

        

        subplot(2,2,subTypes);

        bar(graphTable,'stack')

        title(titleSt);

        hold on;

        

        %plot the propMatSubsetVal1

        subplot(2,2,subTypes+2);

        bar(propMat(:,toPlot'),'stack');

        ylim([0,1]);

        title('(by proportion)')

        if subTypes == 1

            leg = legend(legendText);

            %set(leg,'Location','NorthEastOutside');

            hold on;

        end

    end

    

end



%make outputs for csv printer

header = [header1,header2];

csvData = [graphTable];

output = {header, csvData};