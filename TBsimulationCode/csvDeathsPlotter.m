function csvDeathsPlotter(fullPlottingMatrix, folderName, plotResolution)
%function csvDeathsPlotter(plottingMatrix, plotResolution)
%
%plottingMatrix is a (number of time periods) x 17 matrix for plotting,
%made from TBsimulation with:
%  col 1: nonsmoking natural deaths
%  col 2: smoking natural deaths
%  col 3: nonsmoking untreatedTB deaths
%  col 4: smoking untreatedTB deaths
%  col 5: catI nonsmoking deaths 
%  col 6: catI smoking deaths
%  col 7: catII nonsmoking deaths
%  col 8: catII smoking deaths
%  col 9: catIV nonsmoking deaths
%  col 10: catIV smoking deaths
%  col 11: nonsmoking age-related deaths (reach age 99)
%  col 12: smoking age-related deaths (reach age 99)
%  col13 - 17: all active mdrTB related deaths (see code in TBsimulation_july6 if you
%  need details)
%  col 18: all active sensTB related deaths
%
%folderName is where you want the graphs to print
%plotResolution should be '-r70'


%--------------%%GRAB CHARACTERISTIC COLUMNS%%-------------------%

allSelector = true(1,17);
noneSelector = false(1,17);

numDeathsSelector = noneSelector;
numDeathsSelector(1:12) = true;

smokingSelector = allSelector;
smokingSelector(1:2:17) = false; %turn even columns to false
nonSmokingSelector = ~smokingSelector;

naturalSelector = noneSelector;
naturalSelector(1:2) = true;

untreatedTBSelector = noneSelector;
untreatedTBSelector(3:4) = true;

catIselector = noneSelector;
catIselector(5:6) = true;

catIIselector = noneSelector;
catIIselector(7:8) = true;

catIVselector = noneSelector;
catIVselector(9:10) = true;

oldAgeSelector = noneSelector;
oldAgeSelector(11:12) = true;

mdrSelector = noneSelector;
mdrSelector(13:17) = true;

sensSelector = noneSelector;
sensSelector(18) = true;

%--------------%%MAKE GRAPHS%%-------------------%

for subset = 1:9
    
    if subset == 1
        subset = allSelector;
        titleSuffix = '';
    elseif subset == 2
        subset = (smokingSelector);
        titleSuffix = '_smokingOnly';
    elseif subset == 3
        subset = (nonSmokingSelector);
        titleSuffix = '_nonSmokingOnly';
    elseif subset == 4
        subset = (catIselector);
        titleSuffix = '_CatIOnly';
    elseif subset == 5
        subset = (catIIselector);
        titleSuffix = '_CatIIOnly';
    elseif subset == 6
        subset = (catIVselector);
        titleSuffix = '_CatIVOnly';
    elseif subset == 7
        subset = (sensSelector);
        titleSuffix = '_sensTBonly';
    elseif subset == 8
        subset = (mdrSelector);
        titleSuffix = '_mdrTBonly';
    elseif subset == 9
        subset = (oldAgeSelector);
        titleSuffix = '_oldAgeOnly';
    end

    %%%%%%%%%%%%%health outcomes graph%%%%%%%%%%
    deadTots = sum(fullPlottingMatrix(:,numDeathsSelector),2);
    subsetTots = sum(fullPlottingMatrix(:,subset),2);
    
    toPlot = subsetTots;
    propMat = toPlot ./ deadTots;
    titleSt = strcat('NumberDead', titleSuffix);
    tableTitle = strcat(titleSt, '.png');
    tableHeader = 'number of dead in subset, prop of total dead';
    tableMat = [toPlot,propMat];

    %PLOT
    subplot(2,1,1);
    bar(toPlot)
    title(titleSt);
    hold on;
    subplot(2,1,2);
    bar(propMat);
    ylim([0,1]);
    title('(by proportion of total dead)')
    print('-dpng',plotResolution,[folderName '/' tableTitle]);
    tablePrinter(tableHeader, tableMat, titleSt, folderName);
    close all

end