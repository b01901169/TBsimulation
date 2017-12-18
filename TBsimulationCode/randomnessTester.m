folderNameOut = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\randomTest',

folderName{1} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_fullCatIV_empUptake\2012-01-09_00-02-55',
folderName{2} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_fullCatIV_womUptake\2012-01-09_14-49-22',
folderName{3} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_halfCatIV_empUptake\2012-01-09_08-25-29',
folderName{4} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_halfCatIV_womUptake\2012-01-10_03-45-10',
folderName{5} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_noCatIV_empUptake\2012-01-09_02-08-35',
folderName{6} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\aveTreat_inf0p0020_lat2p55_noCatIV_womUptake\2012-01-09_21-16-56',
folderName{7} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_fullCatIV_empUptake\2012-01-08_19-51-33',
folderName{8} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_fullCatIV_womUptake\2012-01-09_16-57-52',
folderName{9} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_halfCatIV_empUptake\2012-01-09_10-32-37',
folderName{10} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_halfCatIV_womUptake\2012-01-10_05-54-35',
folderName{11} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_noCatIV_empUptake\2012-01-09_04-14-34',
folderName{12} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\bestTreat_inf0p0020_lat2p55_noCatIV_womUptake\2012-01-09_23-26-08',
folderName{13} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\noTreat_inf0p0020_lat2p55_fullCatIV_empUptake\2012-01-08_17-45-42',
folderName{14} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_fullCatIV_empUptake\2012-01-08_21-57-01',
folderName{15} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_fullCatIV_womUptake\2012-01-09_19-07-14',
folderName{16} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_halfCatIV_empUptake\2012-01-09_12-40-00',
folderName{17} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_halfCatIV_womUptake\2012-01-10_08-03-46',
folderName{18} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_noCatIV_empUptake\2012-01-09_06-19-55',
folderName{19} = 'C:\Users\Sze\Dropbox\TBproject\code\outputs\Jan8\worstTreat_inf0p0020_lat2p55_noCatIV_womUptake\2012-01-10_01-35-32',

currentDirectory = pwd,

for i = 1:19
    cd(folderName{i})
    healthOutcomes{i} = dlmread('HealthOutcomes.csv', ',' ,1,0),
    
    %if i >=2
        difference{i} = healthOutcomes{i}
    %end
    
end

currentDirectory = pwd,

zerosMat = zeros( size(difference{1},1),2),
allDifference = [difference{1} , zerosMat, difference{2} , zerosMat, difference{3} ,...
    zerosMat, difference{4} , zerosMat, difference{5} , zerosMat, difference{6} , zerosMat, difference{7} ,...
    zerosMat, difference{8} , zerosMat, difference{9} , zerosMat, difference{10} , zerosMat, difference{11} ,...
    zerosMat, difference{12} , zerosMat, difference{13} , zerosMat, difference{14} , zerosMat, difference{15} ,...
    zerosMat, difference{16} , zerosMat, difference{17} , zerosMat, difference{18}  , zerosMat, difference{19}
    ];

tablePrinter('headerLine', allDifference, 'testing', folderNameOut)