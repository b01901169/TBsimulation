%TBmac outputs plotted
% inputFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\nov20_2014_postLondon_TBmac\paperFigures\means';
% outputFolder = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\nov20_2014_postLondon_TBmac\paperFigures\means\plots';
inputFolder = fullfile(pwd, 'paperFigures', 'means');
outputFolder = fullfile(pwd, 'paperFigures', 'means', 'plots');
mkdir(outputFolder);

%read in file
cd(inputFolder)
% econResults = dlmread('TB_MAC_econ_results_all.csv', ',' ,17,3);
% epiResults = dlmread('TB_MAC_epi_results_all.csv', ',' ,17,3);

econResults = dlmread('mean_econ.csv', ',' ,17,3);
epiResults = dlmread('mean_epi.csv', ',' ,17,3);
daly1 = dlmread('mean_daly1.csv', ',' ,17,3);
daly2 = dlmread('mean_daly2.csv', ',' ,17,3);
daly3 = dlmread('mean_daly3.csv', ',' ,17,4);

%reshape and plot
for colNum = 1:(size(econResults,2))
    selectedEcon = econResults(:,colNum);
    toPlot = reshape(selectedEcon, 20, []);
    plot([2016:2035], toPlot);
    title = strcat('econ_',num2str(colNum));
    print('-dpng','-r70',[outputFolder '/' title]);
    close all
end

for colNum = 1:(size(epiResults,2))
    selectedEcon = epiResults(:,colNum);
    toPlot = reshape(selectedEcon, 20, []);
    plot([2016:2035], toPlot);
    title = strcat('epi',num2str(colNum));
    print('-dpng','-r70',[outputFolder '/' title]);
    close all
end
for colNum = 1:(size(daly1,2))
    selectedEcon = daly1(:,colNum);
    toPlot = reshape(selectedEcon, 20, []);
    plot([2016:2035], toPlot);
    title = strcat('daly1',num2str(colNum));
    print('-dpng','-r70',[outputFolder '/' title]);
    close all
end
for colNum = 1:(size(daly2,2))
    selectedEcon = daly2(:,colNum);
    toPlot = reshape(selectedEcon, 20, []);
    plot([2016:2035], toPlot);
    title = strcat('daly2',num2str(colNum));
    print('-dpng','-r70',[outputFolder '/' title]);
    close all
end
% for colNum = 1:(size(daly3,2))
%     selectedEcon = daly3(:,colNum);
%     toPlot = reshape(selectedEcon, 20, []);
%     plot([2016:2035], toPlot);
%     title = strcat('daly3',num2str(colNum));
%     print('-dpng','-r70',[outputFolder '/' title]);
%     close all
% end