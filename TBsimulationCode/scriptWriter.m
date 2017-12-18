%scriptwriter

%this is the calibrated version:
% TBsimulation_july6_dec9('.', 'b01',180, 480000,  '-r70', 'NA', 1997,0, 2.16, 0.0023, 0.00085, 1.5, 'aveTreat_aveCatIV_empUpta_slowDST_2013')

cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\jan18_vynnicky');

searchArray = {...
[	2.16, 0.0023, 0.00085, 1.5];...
[	0.3	,	0.0026	,	0.0007	,	1.5	];...
[	0.3	,	0.0026	,	0.0006	,	1.5	];...
[	0.3	,	0.0026	,	0.0005	,	1.5	];...
[	0.3	,	0.0027	,	0.0008	,	1.5	];...
[	0.3	,	0.0027	,	0.0007	,	1.5	];...
[	0.3	,	0.0027	,	0.0006	,	2.5	];...
[	0.3	,	0.0027	,	0.0005	,	1.5	];...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
% Default options for calibration 
the_dot = '.'; % '.'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fold_prefix = 0; %remember that this should be a number between 1 and 99 (inclusive) for calibrationMaker script
%%%%%%%%%%%%%%%%%%%%
numYrs = 140;  %numYrs = 180;
numInitPpl = 480000;
graphRes = 70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burnInCom = 'NA';
%%%%%%%%%%%%%%%%%%%%
overWriteYr = 2013;
overWriteYr2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOI_calParam = 2.9;
% activ_calParam = 0.0016;
% uptak_calParam = 0.0031;
% uptakCat2_calParam = 2.4;
%%%%%%%%%%%%%%%%%%%
overWriterFoldName =  'aveTreat_aveCatIV_empUpta_slowDST_2013'; %'base_60perc_private0p5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nIterations = size(searchArray,1);
for iterationNum = 4:nIterations
    fold_prefix = iterationNum;
    FOI_calParam = searchArray{iterationNum}(1);
    activ_calParam = searchArray{iterationNum}(2);
    uptak_calParam = searchArray{iterationNum}(3);
    uptakCat2_calParam = searchArray{iterationNum}(4);
    
    %write the wrapper script
    % Open the file for writing (clear any previous content)
    fileNameStr = sprintf('bar_TBsimulationWrapper%i.m',fold_prefix)
    fid = fopen(fileNameStr,'w')

    % rng('shuffle')
    fprintf(fid,'rng(''shuffle'');\n');

    % Start the function call
    fprintf(fid,'TBsimulation_july6('); % TBsimulation_july6(
    % the dot
    fprintf(fid, '''%s'', ', the_dot); 
    % the output folder name prefix
    fprintf(fid, '''b%0i'', ', fold_prefix); 
    %number of years
    fprintf(fid, '%i, ', numYrs); 
    %num initial people
    fprintf(fid, '%i, ', numInitPpl); 
    %graph res
    fprintf(fid, '''-r%i'', ', graphRes); 
    %burnIn load command
    fprintf(fid, '''%s'', ', burnInCom); 
    %overwrite years
    fprintf(fid, '%i, ', overWriteYr); 
    fprintf(fid, '%i, ', overWriteYr2); 
    %calibration parameters
    fprintf(fid, '%g, ', FOI_calParam); 
    fprintf(fid, '%g, ', activ_calParam); 
    fprintf(fid, '%g, ', uptak_calParam); 
    fprintf(fid, '%g, ', uptakCat2_calParam); 
    %overwriter folder
    fprintf(fid, '''%s'');\n ', overWriterFoldName); 

    %close the script writer
    fclose(fid)
end
