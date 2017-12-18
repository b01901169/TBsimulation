
%to run the end horizon stuff, you have to have the right file in the
%outputs folder.  Go to TBsimParams, find the section on LEbuilder, and run
%that section of the code with the appropriate lines uncommented.  Copy and
%paste the output into the output folder excel files, then run this.
clear
clear all
close all
clc

%%%%%% USER INPUTS  %%%%%%%
stopAndDrop = 0;
makeLEmakerGraphs = 0;
makeCEplane = 1;
makeTwowayGraphs = 0;

takeAverages = 1;

subFolderList = {'base','GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};

%%%%%%%%%%%%%%
% sensitivityStr = 'r28';  %worst GeneX test char
% sensitivityStr = 'r29';  %best GeneX test char
% useAltSensStrFolds = {'base','PPM0p7_0p6_2014'};
% altSensitivityStr = 'p01';
%%%%%%%%%%%%%%
% sensitivityStr = 'r55';  %RNTCP death and default x2 in all interventions 
% useAltSensStrFolds = {'base'};
% altSensitivityStr = 'p01';
%%%%%%%%%%%%%%
% sensitivityStr = 'r56';  %RNTCP death and default x2 in all interventions 
% useAltSensStrFolds = {};
% altSensitivityStr = 'p01';
%%%%%%%%%%%%%%
sensitivityStr = 'r57';  %RNTCP death and default x2 in all interventions 
useAltSensStrFolds = {'base','GeneXallRNTCP_2014','GeneXDST_2014'};
altSensitivityStr = 'p01';
%%%%%%%%%%%%%%
% sensitivityStr = 'r58';  %RNTCP death and default x2 in all interventions 
% useAltSensStrFolds = {'base','GeneXallRNTCP_2014','GeneXDST_2014'};
% altSensitivityStr = 'p01';
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masterFolder = pwd;

fullFolderList = fullfile(masterFolder,subFolderList);
horizonTables = cell(size(subFolderList,1),1);
costAndQalyOutBig = [];

cd(masterFolder)

for i = 1:size(subFolderList,2)
    baseFolder = fullFolderList{i};
    %grab the correct folder(s)
    cd(masterFolder);
    folderNameStartsWith = sensitivityStr;  %this is the first letter of the runs
    
    %check if should use altSensitivtyStr
    useAltSensStrFlag = 0;
    for num = 1:size(useAltSensStrFolds,2)
        if strcmpi(baseFolder,fullfile(masterFolder,useAltSensStrFolds{num}))
            useAltSensStrFlag = 1;
        end
    end
    if useAltSensStrFlag == 1
        folderNameStartsWith = altSensitivityStr;
    end
    
    %for the 10-year-analysis-period
    mustHaveCSV = 'numQalyPpl_postBurnIn.csv';
    %folderNameStr = findMostRecentLegitFolder(baseFolder{i}, folderNameStartsWith, mustHaveCSV)
    folderNameArray = findAllLegitRuns(baseFolder, folderNameStartsWith, mustHaveCSV);  % returns {mostRecent,allLegit}
    %folderNameStr = folderNameArray{1};  %most recent legit run
    sensitivityStrList = folderNameArray{2};   %all legit runs.  Call by allLegitRuns(1) or allLegitRuns(2) etc
    
    %do the same for LEbuilder
    mustHaveCSV = 'numQalyPpl_males.csv';
    folderNameArray = findAllLegitRuns(baseFolder, folderNameStartsWith, mustHaveCSV);  % returns {mostRecent,allLegit}
    %folderNameStr = folderNameArray{1};  %most recent legit run
    LEbuilderSensitivityStr = folderNameArray{2};   %all legit runs.  Call by allLegitRuns(1) or allLegitRuns(2) etc
    
    %reset matracies
    discountedCosts = zeros(1,size(sensitivityStrList,1) );
    discountedQALYs = zeros(1,size(sensitivityStrList,1) );
    totalDisCosts  = zeros(1,size(sensitivityStrList,1) );
    totalDisQalys = zeros(1,size(sensitivityStrList,1) );
    lastState = cell(1,size(sensitivityStrList,1) );
    interpolatedHorizon_costs = zeros(4000,size(LEbuilderSensitivityStr,1));
    interpolatedHorizon_qalys = zeros(4000,size(LEbuilderSensitivityStr,1));
    stDevLEbuilder_costs = zeros(4000,1);
    stDevLEbuilder_qalys = zeros(4000,1);
    
    if stopAndDrop == 1
        size(LEbuilderSensitivityStr,1) = 1;
    end
    
    for LEsensitivityStrNum = 1:size(LEbuilderSensitivityStr,1)
        %GET THE HORIZON TABLE
        if stopAndDrop == 1
            interpolatedHorizon = zeros(4000,3);
        else
            %grab the correct folder
            LEfolderNameStr = LEbuilderSensitivityStr{LEsensitivityStrNum}
            LEbuilderN = getLEbuilderN(LEfolderNameStr);
            %grab horizon table
            [horizonTables,scrapMat, cohortAges, LEbuilderNumCohorts] = qalyMaker2(2014, 2124, LEfolderNameStr, 1, 2014, LEbuilderN);  %first col is LE,2rd col is costs,3nd col is qalys
            
            if makeLEmakerGraphs == 1 && LEsensitivityStrNum == 1 && i == 1
                LEmakerGraphs  %base must be the first subfolder.  Makes debugging graphs
            end
            
            %interpolate the ages, if not all ages were done
            cohorted = reshape(horizonTables, size(cohortAges,1),(LEbuilderNumCohorts/size(cohortAges,1))*3);  %times 3 for LE, qalys, and costs
            interpolated = interp1(cohortAges, cohorted, [0:1:99]', 'linear', 'extrap');
            interpolatedHorizon{1,LEsensitivityStrNum} = reshape(interpolated,(100*(LEbuilderNumCohorts/size(cohortAges,1))),3);
            interpolatedHorizon_costs(:,LEsensitivityStrNum) = interpolatedHorizon{1,LEsensitivityStrNum}(:,2);
            interpolatedHorizon_qalys(:,LEsensitivityStrNum) = interpolatedHorizon{1,LEsensitivityStrNum}(:,3);
        end
        
        %take the average of the LEbuilder horizon costs and qalys for every age-sex-health-trt type
        if LEsensitivityStrNum == size(LEbuilderSensitivityStr,1)
            aveLEbuilder_costs = mean(interpolatedHorizon_costs,2); %should be a 4000 x 1 vector giving ave LEhorizon costs
            aveLEbuilder_qalys = mean(interpolatedHorizon_qalys,2);
            if LEsensitivityStrNum > 1  
                stDevLEbuilder_costs = std(interpolatedHorizon_costs')' / sqrt(size(LEbuilderSensitivityStr,1));  %standard error as a measure of uncertainty of our estimate
                stDevLEbuilder_qalys = std(interpolatedHorizon_qalys')' / sqrt(size(LEbuilderSensitivityStr,1));
            end
        end  %end last LEbuilder run
    end  %end LEsensitivity scenarios loop
    
    for sensitivityStrNum = 1:size(sensitivityStrList,1)
        %GET STOP AND DROP COSTS AND QALYS
        %grab costs and qalys
        cd(masterFolder)
        folderNameStr = sensitivityStrList{sensitivityStrNum}
        [discountedCosts(1,sensitivityStrNum),discountedQALYs(1,sensitivityStrNum)] = qalyMaker2(2014, 2024, folderNameStr, 0,0);
        %grab last state
        cd(folderNameStr)
        lastState{1,sensitivityStrNum} = dlmread(strcat('lastState', '.csv'), ',' ,1,0);
        
        %attach the horizon QALYs and costs to the policy analysis
        lastStateFrac = lastState{1,sensitivityStrNum}/sum(lastState{1,sensitivityStrNum});
        
        %should it be lastStateFrac or lastState{1,sensitivityStrNum}
        totalDisCosts(1,sensitivityStrNum) = discountedCosts(1,sensitivityStrNum) + ( lastStateFrac * aveLEbuilder_costs ); %already discounted when we rolled up
        totalDisQalys(1,sensitivityStrNum) = discountedQALYs(1,sensitivityStrNum) + ( lastStateFrac * aveLEbuilder_qalys );
        
        ave_stDevLEbuilder_costs = (lastStateFrac*stDevLEbuilder_costs);
        ave_stDevLEbuilder_qalys = (lastStateFrac*stDevLEbuilder_qalys);
        
        if sensitivityStrNum == size(sensitivityStrList,1)%grab the means and stdev of the lifetime costs and qalys
            %stop and drop costs and qalys
            aveDiscountedCosts(i,1) = mean(discountedCosts(1,:));
            aveDiscountedQALYs(i,1) = mean(discountedQALYs(1,:));
            stDevDiscountedCosts(i,1) =  std(discountedCosts(1,:))/sqrt(size(sensitivityStrList,2))  ;
            stDevDiscountedQALYs(i,1) =  std(discountedQALYs(1,:))/sqrt(size(sensitivityStrList,2))  ;
            
            %lifetime costs and qalys
            aveLifetimeDiscountedCosts(i,1) = mean(totalDisCosts(1,:));
            aveLifetimeDiscountedQALYs(i,1) = mean(totalDisQalys(1,:));
            stDevLifetimeDiscountedCosts(i,1) = sqrt(  ( stDevDiscountedCosts(i,1)^2  ) + (ave_stDevLEbuilder_costs^2)  );
            stDevLifetimeDiscountedQALYs(i,1) = sqrt(  ( stDevDiscountedQALYs(i,1)^2  ) + (ave_stDevLEbuilder_qalys^2)  );
        end
    end  %end sensitivityStrNum loop
end  %end scenarios loop

% for sensitivityStrNum = 1:size(sensitivityStrList,1)
%     costAndQalyOut = [[1:1:size(subFolderList,2)]',totalDisCosts(:,sensitivityStrNum), totalDisQalys(:,sensitivityStrNum)];
%     costAndQalyOutBig = [costAndQalyOutBig; zeros(1,size(costAndQalyOut,2));costAndQalyOut];
% end


%%%%%%%%%%%%%% print %%%%%%%%%%%%%%%%%%%%%%%
% tableHeader = {'scenarioNum, total costs, total QALYs, all with end horizon.  Rows are, '};
% for j = 1:size(subFolderList,2)
%     tableHeader = strcat(tableHeader(1,1),{','}, subFolderList{j});
% end
% tableHeader = strcat(tableHeader(1,1),{', AND sets are, '});
% for j = 1:size(sensitivityStrList,1)
%     tableHeader = strcat(tableHeader(1,1),{','}, sensitivityStrList{j});
% end
%
% tablePrinter(tableHeader{1}, costAndQalyOutBig, 'costAndQalys', fullfile(masterFolder,'paperFigures'));

%print the average costs and qalys
tableHeader = {'scenarioNum, stopAndDrop aveCosts, stdev, stopAndDrop aveQalys, stdev, lifetime aveCosts, stdev, lifetimeQalys, stdev, Rows are, '};
for j = 1:size(subFolderList,2)
    tableHeader = strcat(tableHeader(1,1),{','}, subFolderList{j});
end
%stop and drop average costs (stdev) and qalys (stdev)
aveCostsAndQalys = [[1:1:size(subFolderList,2)]',aveDiscountedCosts,stDevDiscountedCosts,...
    aveDiscountedQALYs,stDevDiscountedQALYs,aveLifetimeDiscountedCosts,stDevLifetimeDiscountedCosts,...
    aveLifetimeDiscountedQALYs,stDevLifetimeDiscountedQALYs];
cd(masterFolder)
tablePrinter(tableHeader{1}, aveCostsAndQalys, 'aveCostsAndQalys', fullfile(masterFolder,'paperFigures'));

%%%%%%%%%%%%%% make CE plane %%%%%%%%%%%%%%%%%%%%%%%
if makeCEplane == 1
    cd(masterFolder)
    CEplaneMaker
end

%%%%%%%%%%%%%% make twoway graphs if needed %%%%%%%%%%%%%%%%%%%%%%%
if makeTwowayGraphs == 1
    cd(masterFolder)
    twowayPlotMaker
end


