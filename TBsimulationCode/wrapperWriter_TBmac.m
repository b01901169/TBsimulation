%scriptwriter

cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\july8_2014_TBMAC');

overWriterFoldNameList = {...
'base','TBMac1_increaseDots','TBMac1_increaseDots_0p25','TBMac1_increaseDots_0p5','TBMac1_increaseDots_0p75','TBMac1_increaseDots_adv',...
'TBMac2a_improveDotsDefault','TBMac2a_improveDotsDefault_0p25','TBMac2a_improveDotsDefault_0p5','TBMac2a_improveDotsDefault_0p75','TBMac2a_improveDotsDefault_adv',...
'TBMac2b_improveDotsDefault','TBMac2b_improveDotsDefault_0p25','TBMac2b_improveDotsDefault_0p5','TBMac2b_improveDotsDefault_0p75','TBMac2b_improveDotsDefault_adv',...
'TBMac2c_improveDotsDefault','TBMac2c_improveDotsDefault_0p25','TBMac2c_improveDotsDefault_0p5','TBMac2c_improveDotsDefault_0p75','TBMac2c_improveDotsDefault_adv',...
'TBMac3_XpertReplacesSmear','TBMac3_XpertReplacesSmear_0p25','TBMac3_XpertReplacesSmear_0p5','TBMac3_XpertReplacesSmear_0p75','TBMac3_XpertReplacesSmear_adv',...
'TBMac3b_XpertReplacesSmear_SSsensUp','TBMac3b_XpertReplacesSmear_SSsensUp_0p25','TBMac3b_XpertReplacesSmear_SSsensUp_0p5','TBMac3b_XpertReplacesSmear_SSsensUp_0p75','TBMac3b_XpertReplacesSmear_SSsensUp_adv',...
'TBMac4_activeCaseFinding','TBMac4_activeCaseFinding_0p25','TBMac4_activeCaseFinding_0p5','TBMac4_activeCaseFinding_0p75','TBMac4_activeCaseFinding_adv',...
'TBMac5_ACFandLatentTrt_adv',...
'TBMac6_combination','TBMac6_combination_0p25','TBMac6_combination_0p5','TBMac6_combination_0p75','TBMac6_combination_adv'
};
% overWriterFoldNameList = {
% 'TBMac2_improveDotsDefault','TBMac2_improveDotsDefault_0p25','TBMac2_improveDotsDefault_0p5','TBMac2_improveDotsDefault_0p75','TBMac2_improveDotsDefault_adv',...
% 'TBMac2a_improveDotsDefault','TBMac2a_improveDotsDefault_0p25','TBMac2a_improveDotsDefault_0p5','TBMac2a_improveDotsDefault_0p75','TBMac2a_improveDotsDefault_adv',...
% 'TBMac2b_improveDotsDefault','TBMac2b_improveDotsDefault_0p25','TBMac2b_improveDotsDefault_0p5','TBMac2b_improveDotsDefault_0p75','TBMac2b_improveDotsDefault_adv',...
% 'TBMac2c_improveDotsDefault','TBMac2c_improveDotsDefault_0p25','TBMac2c_improveDotsDefault_0p5','TBMac2c_improveDotsDefault_0p75','TBMac2c_improveDotsDefault_adv',...
% };

overWriterFoldNameList = {
'TBMac2_improveDotsDefault','TBMac2_improveDotsDefault_0p25','TBMac2_improveDotsDefault_0p5','TBMac2_improveDotsDefault_0p75','TBMac2_improveDotsDefault_adv',...
};

overWriterFoldNameList = {'TBMac1_increaseDots_adv','TBMac2b_improveDotsDefault_adv','TBMac2c_improveDotsDefault_adv','TBMac2_improveDotsDefault_adv',...
    'TBMac3_XpertReplacesSmear_adv','TBMac3b_XpertReplacesSmear_SSsensUp_adv','TBMac6_combination_adv'};
    

% overWriterFoldNameList = {'base','GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};
% [54,55,69,70]  0
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};  %all geneX scenarios
% [28,29,56,65:68]  24
% overWriterFoldNameList = {'PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'}; % all PPM scenarios
% [30:1:53] 52
% overWriterFoldNameList = {'base','GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};
%[71,72]  use loadStr = 'loadBurnIn';  124

% overWriterFoldNameList = {'PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'}; % all PPM scenarios
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};  %all geneX scenarios
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014','PPM0p7_0p6_2014','PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};
% overWriterFoldNameList = {'GeneXallRNTCP_2014','GeneXDST_2014'};
% overWriterFoldNameList = {'PPM0p7_0p6_GeneXallRNTCP_2014','PPM0p7_0p6_GeneXDST_2014'};

doP01runs = 1;

% loadStr = 'loadBurnInPostTrt';
loadStr = 'loadBurnIn';

barScriptNum = 1;


for folderNum = 1:size(overWriterFoldNameList,2)  %use this folderNumForloop if you want a separate wrapper for every trt arm
    for number = 1:1  %this is the r02 or r03 etc number
        barScriptNum = barScriptNum + 1;  %makes as many wrapper scripts as needed
        
        % Open the file for writing (clear any previous content)
        fileNameStr = sprintf('bar_TBsimulationWrapper%i.m', barScriptNum );
        fid = fopen(fileNameStr,'w');
        
        %initialize
        fold_prefix = strcat('r',sprintf('%02d', number));
        fold_prefixLE = strcat('r',sprintf('%02d', number));
        if doP01runs == 1
            fold_prefix = strcat('p',sprintf('%02d', 1));  %write base
            fold_prefixLE = strcat('p',sprintf('%02d', 1));
        end
        %write the rng('shuffle')
        %fprintf(fid,'rng(''shuffle'');\n');
        sd = floor(29471901*rand);
        fprintf(fid,'rng(%i);\n',sd);
        
        %  for folderNum = 1:size(overWriterFoldNameList,1) %use this folderNumForloop if you want all trt arms in one wrapper
        overWriterFoldName = overWriterFoldNameList{folderNum};
        
        %EXAMPLE  TBsimulation_july6(folderName, logComment,durationYrs, numberPpl, plotResolution, loadBurnInStr, startScenarioYr,startScenarioYr2, latToAct_cal, oldActSlope_cal, oldActIntercept_cal, FOI_cal, aveUptake_cal, cat2uptake_cal, simParamsFolder)
        %EXAMPLE TBsimulation_jan23('.','b01', 160, 200,    '-r70','NA'         , 2014,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3, 'base')
        
        runStr{folderNum} = strcat('TBsimulation_jan23(''.'',''', fold_prefix, ''', 170, 200000, ''-r70'',''', loadStr, ''', 2015,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,''', overWriterFoldName,''');\n');
        fprintf(fid,runStr{folderNum});

        
    end
    %close the script writer
    fclose(fid);
    %end
    
end