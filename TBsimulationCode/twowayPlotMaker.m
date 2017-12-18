%twowayPlotMaker
ICERvalue = [];
sensitivityScenarioCounter = 1;
outputFolder = ('C:\Users\ssuen\Dropbox\TBproject\code\outputs\apr29_2014\paperFigures\sensitivityAnalyses\GeneXcost and DSTcost');

for sensitivityStrNum = 1:size(sensitivityStrList,2)
    sensitivityStr = sensitivityStrList{sensitivityStrNum};
    
    sensitivityScenarioNum = str2num(sensitivityStr(2:3));
    sensitivityScenario;
    xValues(sensitivityStrNum,1) = costs(1, sensitivityScenarioNum);  %these are GeneXcosts
    yValues(sensitivityStrNum,1) = costs(2, sensitivityScenarioNum);  %these are DSTcosts
    
    incrementSize = (size(subFolderList,2)+1);
    ICERvalue(sensitivityScenarioCounter,1) = costAndQalyOutBig(incrementSize*sensitivityScenarioCounter, 6);
    
    sensitivityScenarioCounter = sensitivityScenarioCounter+1;
end


%3D plot
%scatter3(xValues, yValues, ICERvalue);
figure
numPts = 50;
[XI,YI] = meshgrid(linspace(min(xValues),max(xValues),numPts), linspace(min(yValues),max(yValues),numPts));
% now interpolate - find z values for these grid points
ZI = griddata(xValues,yValues,ICERvalue,XI, YI);
mesh(XI,YI,ZI);
xlabel('GeneXCost')
ylabel('DSTcosts')
zlabel('ICERs')
print('-dpng','-r100',[outputFolder '/' 'three_Dplot']);

%2D plots
color{1} = [0,1,0];  %green is go, below WTP
color{2} = [1,0,0];  %red
for WTP = 0:1000:20000
    figure
    for ptNum = 1:size(xValues,1)
        plot(xValues(ptNum), yValues(ptNum),'o','Color',color{1+int8((ICERvalue(ptNum) > WTP))} );
        hold on;
    end
    xlabel('GeneXCost')
    ylabel('DSTcosts')
    title(strcat('Green is Cost Effective with WTP of ',num2str(WTP)));
    print('-dpng','-r100',[outputFolder '/' num2str(WTP)]);
end


mesh(XI,YI,ZI);
xlabel('GeneXCost')
ylabel('DSTcosts')
zlabel('ICERs')
print('-dpng','-r100',[outputFolder '/' 'three_Dplot_flattened']);

