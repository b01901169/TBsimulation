function outMat = recNoDecoder(recNo)
% recNo decoder.  Pass it a column of recNo and it generates the stateMat
% for that time period.  Does not give PageIncr (3), PtrtmntCounter(6) or
% PpastTreatment(9)

Psex = 1; %1 is male, 2 is female
Page = 2;  %age
PageIncr = 3;  %number of months since last bday
Psmoking = 4; %0 is not smoke, 1 is smoke
PurbanRur = 5; %1 is urban, 2 is rural
PtrtmtCounter = 6;  %dummy for counting up months in treatment categories.
Phealth = 7;  % 0=healthy, 1=latentSensTB, 2=latentMDR, 3=activesensTB, 4=activeMDR, 5=unborn, 6=dead
Ptreatment = 8;  %current treatment is 0=none, 1=catI, 2=catII, 3=catIII, 4=catIV
PpastTreatment = 9;  %0=none, 1=ever and noDefault, 2=ever and defaulted

stateMatR(:,Ptreatment) = mod(recNo,6);
recNo = (recNo - stateMatR(:,Ptreatment))/6;
stateMatR(:,Phealth) = mod(recNo,8);
recNo = (recNo - stateMatR(:,Phealth))/8;
stateMatR(:,Psex) = mod(recNo,3);
recNo = (recNo - stateMatR(:,Psex))/3;
stateMatR(:,Page) = mod(recNo,101);
recNo = (recNo - stateMatR(:,Page))/101;
stateMatR(:,Psmoking) = mod(recNo,3);
recNo = (recNo - stateMatR(:,Psmoking))/3;
stateMatR(:,PurbanRur) = mod(recNo,3);

outMat = stateMatR;