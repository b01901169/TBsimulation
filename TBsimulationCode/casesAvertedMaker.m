function casesAvertedMaker(masterFold, comparisonFlag, noTreatRun, comparisonRun, endYear, noteStr)
% incidence comparison generator
% masterFold is a string for the folder in which all the runs live (usually a date.
% Ex:'C:\Users\Sze\Dropbox\TBproject\code\outputs\Feb3\'  )
%
% comparisonFlag is 0 is no comparison, 1 if want to do comparison.
% if comparisonFlag is 1, need the folders to run  (noTreatRun and
% comparisonRun need to be strings to folders.  Ex:
% 'inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake\2012-02-03_01-29-01'  )
%
% noteStr is a string to mark which comparison this was.  Ex:
% 'aveTreatVsNone'


%for endYearNum = [2013, 2038, 2023]
if endYear == 2013;
    endMonthIndex = 216;
    endYearStr = '2013';
elseif endYear == 2038;
    endMonthIndex = 516;
    endYearStr = '2038';
elseif endYear == 2023;
    endMonthIndex = 336;
    endYearStr = '2023';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load incidence output from simulation
currentDirectory = pwd;
comparisonFold = strcat(masterFold, comparisonRun);
cd(comparisonFold);
aveTreatInc = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
aveTreatDeaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
if comparisonFlag == 1
    noTreatFold = strcat(masterFold, noTreatRun);
    cd(noTreatFold);
    noTreatInc = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
    noTreatDeaths = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
end
cd(currentDirectory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load historic population in India
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

TBparams.IndianPop = interp1(IndianPop(:,1), IndianPop(:,2), [1996:1/12:endYear+(11/12)])';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale up the incidences by the scaling factor then multiply by the real pop
aveTreat_numIncInRealPop(:,1:4) =( aveTreatInc(1:endMonthIndex,1:4) ./ repmat(aveTreatInc(1:endMonthIndex,5),1,4) ) .* repmat(TBparams.IndianPop,1,4);
aveTreat_numDeathsInRealPop(:,1:2) =( aveTreatDeaths(1:endMonthIndex,1:2) ./ repmat(aveTreatDeaths(1:endMonthIndex,3),1,2) ) .* repmat(TBparams.IndianPop,1,2);

%find the difference in the incidences

if comparisonFlag == 1
    noTreat_numIncInRealPop(:,1:4) = ( noTreatInc(1:endMonthIndex,1:4) ./ repmat(noTreatInc(1:endMonthIndex,5),1,4) ) .* repmat(TBparams.IndianPop,1,4);
    noTreat_numDeathsInRealPop(:,1:2) =( noTreatDeaths(1:endMonthIndex,1:2) ./ repmat(noTreatDeaths(1:endMonthIndex,3),1,2) ) .* repmat(TBparams.IndianPop,1,2);
    casesAvertedVec = sum( noTreat_numIncInRealPop - aveTreat_numIncInRealPop ) ;
    deathsAvertedVec = sum( noTreat_numDeathsInRealPop - aveTreat_numDeathsInRealPop ) ;
    
    header = 'averted latent DS infection cases, averted MDR infection cases,averted act DS disease cases, averted active MDR disease cases';
    csvFileName = strcat('casesAverted_from1996_', endYearStr,'_',noteStr);
    tablePrinter(header, casesAvertedVec, csvFileName, comparisonFold);
%     disp(csvFileName)
%     disp(header);
%     casesAvertedVec
    
    header = 'averted DS deaths, averted MDR deaths';
    csvFileName = strcat('deathsAverted_from1996_', endYearStr,'_',noteStr);
    tablePrinter(header,deathsAvertedVec, csvFileName, comparisonFold);
%     disp(csvFileName)
%     disp(header);
%     casesAvertedVec
else
    header = 'total incident latent DS infection cases, total incident MDR infection cases,total incident act DS disease cases, total incident active MDR disease cases';
    csvFileName = strcat('casesTotal_from1996_', endYearStr,'_',noteStr);
    tablePrinter(header,  sum(aveTreat_numIncInRealPop,1), csvFileName, comparisonFold);
%     disp(csvFileName)
%     disp(header);
%     sum(aveTreat_numIncInRealPop,1)
    
    header = 'total DS deaths, total MDR deaths';
    csvFileName = strcat('deathsTotal_from1996_', endYearStr,'_',noteStr);
    tablePrinter(header, sum(aveTreat_numDeathsInRealPop,1), csvFileName, comparisonFold);
%     disp(csvFileName)
%     disp(header);
%     sum(aveTreat_numDeathsInRealPop,1)
    
    
end
clearvars aveTreat_numIncInRealPop noTreat_numIncInRealPop casesAvertedVec aveTreat_numDeathsInRealPop noTreat_numDeathsInRealPop deathsAvertedVec
%end  %end endYear loop



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     cd(comparisonFold);
%     deathsPostBurnIn = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
%     diseaseCasesMat = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
%     popDeaths = zeros(size(deathsPostBurnIn,2),2);
%     popCases = zeros(size(deathsPostBurnIn,2), 4);
%     for timePeriod = 1: size(deathsPostBurnIn,2);  %if 1 is 1996, month 193 is jan2012.  492 is dec2036.
%         popMonth = timePeriod;
%         %get deaths
%         popDeaths(timePeriod,:) =  TBparams.IndianPop(popMonth)*( deathsPostBurnIn(timePeriod,1:2) ./ repmat(deathsPostBurnIn(timePeriod,3), 1, 2) );
%
%         %get cases
%         popCases(timePeriod,:) =  TBparams.IndianPop(popMonth)*(  diseaseCasesMat(timePeriod,1:4,:) ./ repmat(diseaseCasesMat(timePeriod,5), 1, 4));
%     end
%
%     cd(noTreatFold);
%     deathsPostBurnIn = dlmread('TBdeaths_postBurnIn.csv', ',' ,1,0);
%     diseaseCasesMat = dlmread('DSmdrIncidence_postBurnIn.csv', ',' ,1,0);
%     popDeathsNoTreat = zeros(size(deathsPostBurnIn,2),2);
%     popCasesNoTreat = zeros(size(deathsPostBurnIn,2), 4);
%     for timePeriod = 1: size(deathsPostBurnIn,2);  %if 1 is 1996, month 193 is jan2012.  492 is dec2036.
%         popMonth = timePeriod;
%         %get deaths
%         popDeathsNoTreat(timePeriod,:) =  TBparams.IndianPop(popMonth)*( deathsPostBurnIn(timePeriod,1:2) ./ repmat(deathsPostBurnIn(timePeriod,3), 1, 2) );
%         %get cases
%         popCasesNoTreat(timePeriod,:) =  TBparams.IndianPop(popMonth)*(  diseaseCasesMat(timePeriod,1:4,:) ./ repmat(diseaseCasesMat(timePeriod,5), 1, 4));
%     end

end

%     currentDirectory = pwd;
%     cd(comparisonFold);
%     aveTreatInc2012_2036 = dlmread('TotalCases_scaledToIndianPop.csv', ',' ,1,0);
%     aveTreatDeaths2012_2036 = dlmread('TotalDeaths_scaledToIndianPop.csv', ',' ,1,0);
%     aveTreatQALY2012_2036 = dlmread('TotalQALYs_scaledToIndianPop.csv', ',' ,1,0);
%     cd(noTreatFold);
%     noTreatInc2012_2036 = dlmread('TotalCases_scaledToIndianPop.csv', ',' ,1,0);
%     noTreatDeaths2012_2036 = dlmread('TotalDeaths_scaledToIndianPop.csv', ',' ,1,0);
%     noTreatQALY2012_2036 = dlmread('TotalQALYs_scaledToIndianPop.csv', ',' ,1,0);
%     cd(currentDirectory);
%
%     casesAverted = sum(noTreatInc2012_2036 )-sum(aveTreatInc2012_2036);
%     deathsAverted = sum(noTreatDeaths2012_2036)-sum(aveTreatDeaths2012_2036);
%     QALYSgained = sum(aveTreatQALY2012_2036)-sum(noTreatQALY2012_2036);
%
%     headerArray{1} = 'averted latent DS infection cases, averted MDR infection cases,averted act DS disease cases, averted active MDR disease cases';
%     headerArray{2} = 'avertedDSdeaths, avertedMDRdeaths';
%     headerArray{3} = 'totalQALYsGained';
%     csvFileName = {'casesAverted_from2012_2036', 'deathsAverted_from2012_2036', 'QALYSgained_from2012_2036'};
%     csvFile = {casesAverted,deathsAverted, QALYSgained };
%     for i = 1:3
%         tablePrinter(headerArray{i}, csvFile{i}, csvFileName{i}, comparisonFold);
%     end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % %    display
% %     fprintf('latDS cases averted is %d\n',casesAvertedVec(1));
% %     fprintf('latMDR cases averted is %d\n',casesAvertedVec(2));
% %     fprintf('actDS cases averted is %d\n',casesAvertedVec(3));
% %     fprintf('actMDR cases averted is %d\n',casesAvertedVec(4));
% %





