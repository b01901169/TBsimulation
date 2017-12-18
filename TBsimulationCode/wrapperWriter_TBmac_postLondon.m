%scriptwriter

%cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\dec13_2014_postLondon_TBmac');
cd('D:\Dropbox\TBmac\June15p5');

TBmac_folders;

% loadStr = 'loadBurnInPostTrt';
loadStr = 'loadBurnIn';
barScriptNum = 0;
doP01runs = 0;
simulationPop = 200000;
% simulationPop = 2000;  %DEBUGGING

for folderNum = 1:size(overWriterFoldNameList,2)  %use this folderNumForloop if you want a separate wrapper for every trt arm
    for number = 3:11  %this is the r02 or r03 etc number
        barScriptNum = barScriptNum + 1;  %makes as many wrapper scripts as needed
        
        % Open the file for writing (clear any previous content)
        fileNameStr = sprintf('bar_TBsimulationWrapper%i.m', barScriptNum );
        fid = fopen(fileNameStr,'w');
        
        %initialize
        fold_prefix = strcat('p',sprintf('%02d', number));
        fold_prefixLE = strcat('p',sprintf('%02d', number));
        if doP01runs == 1
            fold_prefix = strcat('p',sprintf('%02d', 1));  %write base
            fold_prefixLE = strcat('p',sprintf('%02d', 1));
        end
        %write the rng('shuffle')
        %fprintf(fid,'rng(''shuffle'');\n');
        sd = floor(29471901*rand);
        fprintf(fid,'rng(%i);\n',sd);
        fprintf(fid,'cd(''/home/ssuen/TBmac/'');\n');

        %  for folderNum = 1:size(overWriterFoldNameList,1) %use this folderNumForloop if you want all trt arms in one wrapper
        overWriterFoldName = overWriterFoldNameList{folderNum};
        
        %EXAMPLE  TBsimulation_july6(folderName, logComment,durationYrs, numberPpl, plotResolution, loadBurnInStr, startScenarioYr,startScenarioYr2, latToAct_cal, oldActSlope_cal, oldActIntercept_cal, FOI_cal, aveUptake_cal, cat2uptake_cal, simParamsFolder)
        %EXAMPLE TBsimulation_jan23('.','b01', 160, 200,    '-r70','NA'         , 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
        
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefix, ''', 170,', num2str(simulationPop),', ''-r70'',''', loadStr, ''', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum}{1});
        
    end
    %close the script writer
    fclose(fid);
    %end
    
end

fclose('all');