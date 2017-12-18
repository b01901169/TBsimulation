%CEplaneMaker make the CE plane.  Called from LEmaker
%C:\Users\ssuen\Dropbox\TBproject\code\outputs

if takeAverages == 1
    costsAndQalys = aveCostsAndQalys(:,2:end);
    costsCol = 1;
    stdDevCost = 2;
    qalyCol = 3;
    stdDevQaly = 4;

    repeatTimes = 2;
    repeatTimesStr = {'stopAndDrop','Lifetime'};
    subFolderList = {subFolderList{:},subFolderList{:}};
    xMinAdj = [2,10];
else
    costsAndQalys = [costAndQalyOutBig(2:end,2:3), zeros(size(aveCostsAndQalys,1)), zeros(size(aveCostsAndQalys,1))];  %from LEmaker
    costsCol = 1;
    qalyCol = 2;
    stdDevCost = 3;
    stdDevQaly = 4;
    
    repeatTimes = 1;
    repeatTimesStr = 'noAveraging';
    xMinAdj = [2,10];

end

    for i = 1:size(subFolderList,2)
        subFolderListSpaces(1,i) = strrep(subFolderList(1,i), '_', ' ');
    end


for repTime = 1:repeatTimes
    repTime
    if takeAverages == 1 && repTime == 2
        costsCol = 5;
        stdDevCost = 6;
        qalyCol = 7;
        stdDevQaly = 8;
    end
    
    % order by total costs
    [scrap, sortedIndex] = sort(costsAndQalys(:,costsCol));
    sortedCostsQalys = costsAndQalys(sortedIndex,:);
    
    % plot total qalys by total costs on CE plane
    figure
    labels = subFolderListSpaces(:,sortedIndex);  %labels correspond to their order
    cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs')
    errorbarxy(sortedCostsQalys(:,costsCol),  sortedCostsQalys(:,qalyCol), sortedCostsQalys(:,stdDevCost)*1.96,sortedCostsQalys(:,stdDevCost)*1.96, sortedCostsQalys(:,stdDevQaly)*1.96,sortedCostsQalys(:,stdDevQaly)*1.96)
    cd(masterFolder)
    %errorbar(sortedCostsQalys(:,costsCol),  sortedCostsQalys(:,qalyCol),sortedCostsQalys(:,stdDevQaly)*1.96,sortedCostsQalys(:,stdDevQaly)*1.96,'rx');
    %plot(sortedCostsQalys(:,costsCol), sortedCostsQalys(:,qalyCol), 'rx')
    hold on;
    plot(sortedCostsQalys(:,costsCol), sortedCostsQalys(:,qalyCol))
    hold on;
    text(sortedCostsQalys(:,costsCol), sortedCostsQalys(:,qalyCol), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
    xlabel('Total Costs')
    ylabel('Total QALYs')
    title(strcat('CE Plane: All Treatment Arms, ',repeatTimesStr{repTime}) )
    %xlim([min(sortedCostsQalys(:,costsCol))-xMinAdj(repTime) max(sortedCostsQalys(:,costsCol))+5]);
    %ylim([min(sortedCostsQalys(:,qalyCol))-1 max(sortedCostsQalys(:,qalyCol))+1]);
    titleStr = strcat('CE_Plane_all_Treatment_Arms',repeatTimesStr{repTime}) ;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 7])
    print('-dpng','-r300',[fullfile(masterFolder,'paperFigures') '/' titleStr]);
    
    
    disp('removing dominated')
    % removed dominated treatment arms
    delQalys = sortedCostsQalys(:,qalyCol)-[0;sortedCostsQalys(1:end-1,qalyCol)];  %Qalys
    undominatedCostQalys = sortedCostsQalys(delQalys >= 0,:);
    
    % plot undominated total qalys by total costs on CE plane
    figure
    newlabels = labels(:,delQalys >= 0);  %' # labels correspond to their order
    cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs')
    errorbarxy(undominatedCostQalys(:,costsCol),  undominatedCostQalys(:,qalyCol), undominatedCostQalys(:,stdDevCost)*1.96,undominatedCostQalys(:,stdDevCost)*1.96, undominatedCostQalys(:,stdDevQaly)*1.96,undominatedCostQalys(:,stdDevQaly)*1.96)
    cd(masterFolder)    
    %errorbar(undominatedCostQalys(:,costsCol),  undominatedCostQalys(:,qalyCol),undominatedCostQalys(:,stdDevQaly)*1.96,undominatedCostQalys(:,stdDevQaly)*1.96,'rx');
    % plot(undominatedCostQalys(:,costsCol), undominatedCostQalys(:,qalyCol), 'rx')  %points
    hold on;
    plot(undominatedCostQalys(:,costsCol), undominatedCostQalys(:,qalyCol))  %line
    hold on;
    text(undominatedCostQalys(:,costsCol), undominatedCostQalys(:,qalyCol), newlabels, 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
    xlabel('Total Costs')
    ylabel('Total QALYs')
    title(strcat('CE Plane: Undominated Treatment Arms, ',repeatTimesStr{repTime}))
    %xlim([min(undominatedCostQalys(:,costsCol))-xMinAdj(repTime) max(undominatedCostQalys(:,costsCol))+5]);
    %ylim([min(undominatedCostQalys(:,qalyCol))-1 max(undominatedCostQalys(:,qalyCol))+1]);
    titleStr = strcat('CE_Plane_Un_strictly_dominated_Treatment_Arms',repeatTimesStr{repTime}) ;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 7])
    print('-dpng','-r300',[fullfile(masterFolder,'paperFigures') '/' titleStr]);
    
    % make ICERS
    deltaCostQalys = undominatedCostQalys(:,[costsCol,qalyCol]) - [zeros(1,2);undominatedCostQalys(1:end-1,[costsCol,qalyCol])];
    withIcers = [undominatedCostQalys(:,[costsCol,qalyCol]) , deltaCostQalys, deltaCostQalys(:,1)./deltaCostQalys(:,2)];
    
    fullstr = '';
    for ind = 1:size(newlabels,2)
        fullstr = strcat(fullstr, ',',newlabels{1,ind});
    end
    tableHeader = strcat('Total costs, total QALYs, incre costs, incre QALYS, ICERs, undominated strategies only, rows are: ', fullstr);
    tableName = strcat('unDominatedICERs',repeatTimesStr{repTime});
    tablePrinter(tableHeader, withIcers, tableName, fullfile(masterFolder,'paperFigures'));
end

