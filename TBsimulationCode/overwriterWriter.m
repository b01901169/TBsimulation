%overwriter writer

clear
clear all
clc
cd('D:\Dropbox\TBmac\June15p5')

%get folder list
TBmac_folders;

%initialize folders
masterFolder = pwd;
sourceFolder = 'base\';

for scenarioNum = 1:size(scenarioList,1)
    scenarioNameStr = scenarioList{scenarioNum,1}{1};
    multiplierVal = scenarioList{scenarioNum,1}{2};
    scenarioNameRow = {scenarioList{scenarioNum,1}{3:end}};
    cd(masterFolder);
    
    %make the new overwriter directoy
    mkdir(scenarioNameStr{1});
    
    %copy the overwriter and post2011 overwriter
    source = fullfile(masterFolder, sourceFolder, 'TBsimParamsOverwriter*');
    destination = fullfile(masterFolder,scenarioNameStr{1});
    status  = copyfile(source,destination);
    
    %change the post2011 overwriter to turn on the right scenarios
    cd(destination)
    
    %get the text
    fileID = fopen('TBsimParamsOverwriter_post2011.m', 'r');
    modifiedStr = fscanf(fileID,'%c');
    fclose(fileID);
    
    %change each one of the modifiers from 0 to 1
    for scenarioNameVal = 1:size(scenarioNameRow,2)
        scenarioName = scenarioNameRow{scenarioNameVal};
        %check that modification happened
        if strfind(modifiedStr, strcat(scenarioName{1}, ' = 0')) == 0
            disp(scenarioName{1});
            error('did not find string');
        else  %change it
            modifiedStr = strrep(modifiedStr, strcat(scenarioName{1}, ' = 0'), strcat(scenarioName{1}, ' = 1'));
        end
    end
    
    %check that change the multiplier
    if strfind(modifiedStr, 'TBparams.TBmacMultiplier = 1') == 0
        disp('TBparams.TBmacMultiplier = 1');
        error('did not find string');
    else  %change it
        modifiedStr = strrep(modifiedStr, 'TBparams.TBmacMultiplier = 1', strcat('TBparams.TBmacMultiplier = ', num2str(multiplierVal{1})));
    end
        
    %write
    fileID = fopen('TBsimParamsOverwriter_post2011.m', 'w');
    fprintf(fileID,'%s',modifiedStr);
    fclose(fileID);
    
    
end

fclose('all');