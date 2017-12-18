%run this script, then upload the new .scripts from the zzzz_barleyScripts
%folder to your afs space.  Copy and paste the commands in the
%aaa_convertToUnix.txt file to your SSH command and press enter.  Copy and
%paste the commands from qsub.txt into your SSH command and press enter.

unixFolder = 'aug2_flux_'; cornFolderName = 'aug2_flux';
unixFolderShort = 'aug2_flux_';
totNumScriptsFrom = 1;  %1
totNumScriptsTo = 24;

%makes the commands to qsub
cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs')
fileNameStr = 'qsub.txt';
fid = fopen(fileNameStr,'wt');

for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    statement = char(strcat('qsub',{' '},unixFolder, numStr,'.script\n'));
    fprintf(fid,statement);
end
fprintf(fid,'qstat');
fclose(fid);

%makes the .scripts
cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\zzzz_barleyScripts')
for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    fileNameStr = strcat(unixFolderShort, numStr,'.script');
    fid = fopen(fileNameStr,'wt');
    
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,char(strcat('#$ -N', {' '},unixFolder,numStr, '\n')));
    fprintf(fid,'#$ -m be\n');
    fprintf(fid,'#$ -M ssuen@stanford.edu\n');
    fprintf(fid,char(strcat('cd /afs/ir.stanford.edu/users/s/s/ssuen/Matlab/',cornFolderName,'/\n')));
    fprintf(fid,char(strcat('/afs/ir.stanford.edu/software/matlab-2011b/bin/matlab -nodesktop < bar_TBsimulationWrapper',numStr ,'.m')));
    
end
fclose(fid);

%makes the commands to convert the .scripts to unix
cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\zzzz_barleyScripts')
fileNameStr = 'aaa_convertToUnix.txt';
fid = fopen(fileNameStr,'wt');

for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    fprintf(fid,char(strcat('dos2unix',{' '}, unixFolder,numStr ,'.script\n')));
    
end
fclose(fid);
