%LEmakerGraphs
%script used in LEmaker

dir = pwd;
cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs')
pwd
outputFolder = fullfile(masterFolder,'paperFigures');
leStates_raw = dlmread('LEbuilderStates.csv', ',' ,1,1);
cd(dir)
leStates = leStates_raw(:,[1,2,7,8,10,14]);
leStatePlotData = [leStates,horizonTables];

Psex = 1;
Page = 2;
Phealth = 3;
Ptreatment = 4;
PtimePerInfect = 5;
PactTime = 6;
Ple = 7;
Pcosts = 8;
Pqalys = 9;

possibleHealthStates = [0:4];
possibleTrtStates = [0,1,2,4];
sexStr = {'Male ','Female '};
titleStr = {'Uninfected','Latent DS',' Latent MDR','Active DS','Active MDR'};
outputVarStr = {'LE','Costs','QALYs'};
ylimMat = [0,100; 0, 1400 ; 0, 14; ];  %for no discounting
%ylimMat = [0,100; 0, 2*10^13 ; 0, 2*10^10; ];
males = (leStatePlotData(:,Psex) == 1);
females = (leStatePlotData(:,Psex) == 2);
atBirth = find(leStatePlotData(males,Page) == 0);  %index of age=0 cohorts

sexSeparated = {leStatePlotData(males,:),leStatePlotData(females,:)};

outputVar = Ple;
for outputVar = [Ple, Pcosts, Pqalys]
    figure
    counter = 1;
    for sexVal = 1:2
        if outputVar == Ple
            %%%%%%%%%% LE at birth tables %%%%%%%%%%
            outTable = reshape(sexSeparated{sexVal}(atBirth,outputVar),4,5)';
            %format with the values and output
            outTable = [possibleHealthStates',outTable];            
            fprintf('Life expectancy at birth for sex val %i', sexVal)
            outTable = [[0,possibleTrtStates];outTable]
            LEtables{sexVal} = outTable;
        end
        
        %%%%%%%%%% LE graphs %%%%%%%%%%
        graphData = sexSeparated{sexVal};
        
        longGraphData = reshape(graphData(:,outputVar),size(cohortAges,1),20);  %now age increments down each col, and we only have the outputVar
        longGraphData = [longGraphData;zeros(1,size(longGraphData,2))];
        
        for hlthIndexNum = [1:4:17]
            subplot(2,5,counter)
            horizonOut{outputVar, counter} = longGraphData(:,hlthIndexNum:hlthIndexNum+3);  %this is used in the GeneXpertTiming proj
            plot([cohortAges;100],longGraphData(:,hlthIndexNum:hlthIndexNum+3));
            if counter == 1 && outputVar == 7
                legend('No trtmt','Cat I','Cat II','Cat IV');
            end
            title(strcat(sexStr{sexVal},titleStr{ceil(hlthIndexNum/4)}));
           % ylim(ylimMat(outputVar-6,:));
            xlabel('Age')
            plotTitle = outputVarStr{outputVar-6};
            %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
            print('-dpng','-r100',[outputFolder '/' plotTitle]);
            
            counter = counter + 1;
        end
    end
end



