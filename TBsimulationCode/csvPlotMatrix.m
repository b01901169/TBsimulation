function csvPlotMatrix(recNo, folderName)
%function csvPlotMatrix(recNo)
%Makes a csv file with counts of people in different health, treatment,
%sex, smoking, and age bins over all time
%
%RECNO is the compressed stateMatrix with all the information in it.
%csvFile has timePeriod on rows and columns are:
% (male, smoking, age0_20, healthy), (male, smoking, age0_20,latentSens)...
% (male, smoking, age21_60, healthy),(male, smoking, age21_60, latentSens)...
% (male, smoking, age61_100, healthy),...
% (male, nonsmoking, age0_20, healthy),...
% ...(health outcomes change, ages change)
% (female, smoking, age0_20, healthy),...
% ...(health outcomes, ages, and smoking changes, in order given above)
%  There should be 168 columns.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go to the output directory
currentDirectory = pwd;
cd(folderName);

%initialize vars
totPeriods = size(recNo,2);
numCodes = max(max(recNo));

%make the definition of the categories
codeNumMat = zeros(numCodes,1);
%make matrix of code values
for i = 1: numCodes
    codeNumMat(i,1) = i;
end

%--------------%%MAKE HEALTH AND TREATMENT VALUES%%-------------------%

%always need a healthMat to read the health codes
healthMat = int8(floor( (  mod(codeNumMat,48) )/6 ));

%decode the treatment codes
trtmtVal = int8(mod(codeNumMat,6));

%make indicies easier to read;
numHealthOutcomes = 14;
healthy = 1; latentSens = 2; latentMdr = 3; activeSens =4; activeMdr = 5; notYetBorn = 6; dead = 7;
catI_sens = 8; catII_sens = 9; catIII = 10; catIV_sens = 11; catI_MDR = 12; catII_MDR = 13; catIV_MDR = 14;

% Generate healthTreatVal
healthTreatVal = false(numCodes, numHealthOutcomes);
healthTreatVal(:, healthy) = (healthMat == 0 & trtmtVal == 0 );
healthTreatVal(:, latentSens) = (healthMat == 1 & trtmtVal == 0 );
healthTreatVal(:, latentMdr) = (healthMat == 2 & trtmtVal == 0 );
healthTreatVal(:, activeSens) = (healthMat ==3 & trtmtVal == 0 );
healthTreatVal(:, activeMdr) = (healthMat == 4 & trtmtVal == 0 );
healthTreatVal(:, notYetBorn) = (healthMat == 5);
healthTreatVal(:, dead) = (healthMat == 6);
healthTreatVal(:, catI_sens) = (trtmtVal == 1 & healthMat == 3 );
healthTreatVal(:, catII_sens) = (trtmtVal == 2 & healthMat == 3 );
healthTreatVal(:, catIII) = (trtmtVal == 3 & healthMat == 3);
healthTreatVal(:, catIV_sens) = (trtmtVal == 4 & healthMat == 3 );
healthTreatVal(:, catI_MDR) = (trtmtVal == 1 & healthMat == 4 );
healthTreatVal(:, catII_MDR) = (trtmtVal == 2 & healthMat == 4 );
healthTreatVal(:, catIV_MDR) = (trtmtVal == 4 & healthMat == 4 );

%--------------%%MAKE SEX, SMOKING, AGE VALUES%%-------------------%

%make subset matracies
sexMat = int8((1/48)*( (mod(codeNumMat,144)) - 6*double(healthMat) - double(trtmtVal))); %sex
smokingMat = int8(floor(mod(codeNumMat,43632)/14544));  %smoking
ageMat =  int8( floor( mod(codeNumMat,14544)/144 ) );   %age

%make types easier to read
male = 1; female = 2; 
smoking = 1; nonSmoking = 2; 
ageBinMat_1 = 1; ageBinMat_2 = 2; ageBinMat_3 = 3;

%turn sex, smoking, age into logicals that tell me which bin the recNo codes are in
sexLogicalMat = false(numCodes,2);
sexLogicalMat(:,male) = (sexMat == 1); 
sexLogicalMat(:,female) = (sexMat == 2);  

smokingLogicalMat = false(numCodes,2);
smokingLogicalMat(:,smoking) = (smokingMat == 1); 
smokingLogicalMat(:,nonSmoking) = (smokingMat == 0); 

ageLogicalMat = false(numCodes,3);
ageLogicalMat(:,ageBinMat_1) = (ageMat <= 20); %1 if age 0_20
ageLogicalMat(:,ageBinMat_2)  = ( 21 < ageMat & ageMat <= 60); %1 if age 21_60
ageLogicalMat(:,ageBinMat_3) = ( 61 < ageMat & ageMat <= 110); %1 if age 61_100

%--------------%%MAKE INDICATOR MATRIX WITH EACH DEMOG TYPE%%-------------------%

%subsetMat gives a 0 or 1 for whether recNo code is in demographic subset
column = 1;
subsetMat = false(numCodes,(  size(sexLogicalMat,2)*size(smokingLogicalMat,2)*size(ageLogicalMat,2)*numHealthOutcomes   ));
for sexType = 1:size(sexLogicalMat,2);
    for smokingType = 1:size(smokingLogicalMat,2);
        for ageType = 1:size(ageLogicalMat,2);
            for healthTreatType = 1 : numHealthOutcomes;
                subsetMat(:,column) = (sexLogicalMat(:,sexType) & smokingLogicalMat(:,smokingType) & ageLogicalMat(:,ageType) & healthTreatVal(:,healthTreatType) );
                column = column + 1;  %increment a column for each demographic type
            end
        end
    end
end

%--------------%%COUNT NUMBER OF PEOPLE IN EACH TYPE%%-------------------%
csvMat = zeros(totPeriods, size(subsetMat,2), 'int32');

%TIME LOOP
for timePeriod = 1: totPeriods
    %reset vars
    codeFreqMat = zeros(numCodes,1);
    
    % sum over recNo codes to make a table of frequency of each recNo code
    for person = 1 : size(recNo,1)
        code = recNo(person,timePeriod);
        codeFreqMat(code) = codeFreqMat(code) + 1;
    end
    
    for demogType = 1:size(subsetMat,2)
        csvMat(timePeriod, demogType) = sum(codeFreqMat(subsetMat(:,demogType)), 1);
    end
    
end  %end timePeriod loop

%--------------%%OUTPUT A CSV FILE%%-------------------%

%write this to a csvFile
dlmwrite('csvMatForPlotting.txt', csvMat);

%--------------%%CHANGE DIR%%-------------------%

%get back out to the original directory
cd(currentDirectory);