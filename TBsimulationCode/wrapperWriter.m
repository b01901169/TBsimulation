%scriptwriter

cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\july7_2014');

overWriterFoldNameList = {'base','GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};

%overWriterFoldNameList = {'PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'}; % all PPM scenarios
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};  %all geneX scenarios
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};

% loadStr = 'loadBurnInPostTrt';
loadStr = 'loadBurnIn';

barScriptNum = 112;


for folderNum = 1:size(overWriterFoldNameList,2)  %use this folderNumForloop if you want a separate wrapper for every trt arm
    for number = 69:70 %2:26  %this is the r02 or r03 etc number
        barScriptNum = barScriptNum + 1;  %makes as many wrapper scripts as needed
        
        % Open the file for writing (clear any previous content)
        fileNameStr = sprintf('bar_TBsimulationWrapper%i.m', barScriptNum );
        fid = fopen(fileNameStr,'w');
        
        %initialize
        fold_prefix = strcat('r',sprintf('%02d', number));
        fold_prefixLE = strcat('r',sprintf('%02d', number));
        
        %write the rng('shuffle')
        %fprintf(fid,'rng(''shuffle'');\n');
        sd = floor(29471901*rand);
        fprintf(fid,'rng(%i);\n',sd);
        
        %  for folderNum = 1:size(overWriterFoldNameList,1) %use this folderNumForloop if you want all trt arms in one wrapper
        overWriterFoldName = overWriterFoldNameList{folderNum};
        
        %EXAMPLE  TBsimulation_july6(folderName, logComment,durationYrs, numberPpl, plotResolution, loadBurnInStr, startScenarioYr,startScenarioYr2, latToAct_cal, oldActSlope_cal, oldActIntercept_cal, FOI_cal, aveUptake_cal, cat2uptake_cal, simParamsFolder)
        %EXAMPLE TBsimulation_jan23('.','b01', 160, 200,    '-r70','NA'         , 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
        
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefix, ''', 160, 200000, ''-r70'',''', loadStr, ''', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        %five LEbuilder runs
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefix, ''', 160, 200000, ''-r70'',''', loadStr, ''',2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        %five LEbuilder runs
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefixLE, ''', 160, 100000, ''-r70'', ''LEbuilder'', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});
        
    end
    %close the script writer
    fclose(fid);
    %end
    
end