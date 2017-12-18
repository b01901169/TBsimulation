%run this script, then upload the new .scripts from the zzzz_barleyScripts
%folder to your afs space.  Copy and paste the commands in the
%aaa_convertToUnix.txt file to your SSH command and press enter.  Copy and
%paste the commands from qsub.txt into your SSH command and press enter.

unixShellName = 'TBmacW'; sherlockFolderName = 'TBmac';
totNumScriptsFrom = 1;  %1
totNumScriptsTo = 108;

outFoldName = 'D:\Dropbox\TBmac\zzzz_sherlockScripts';

%makes the commands to qsub
% cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\zzzz_sherlockScripts')
cd(outFoldName)
fileNameStr = 'aaa_sbatch.txt';
fid = fopen(fileNameStr,'wt');

for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    statement = char(strcat('sbatch',{' '},unixShellName, numStr,'.script\n'));
    fprintf(fid,statement);
end
fprintf(fid,'sacct');
fclose(fid);

% %makes the .scripts
% #!/bin/bash
% #SBATCH --job-name=TBmac2_job
% #SBATCH --output=TBmac2_job.out
% #SBATCH --error=TBmac2_job.err
% #SBATCH --time=2:00:00
% #SBATCH --qos=normal
% #SBATCH --nodes=1
% #SBATCH --mem=10000
% #SBATCH --mail-type=ALL
% #SBATCH  --mail-user=ssuen@stanford.edu
% #################
% module load matlab
% matlab -nosplash -nodesktop -nodisplay -singleCompThread < $HOME/TBmac/bar_TBsi$

cd(outFoldName)
for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    fileNameStr = strcat(unixShellName, numStr,'.script');
    fid = fopen(fileNameStr,'wt');
    
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,char(strcat('#SBATCH --job-name=',unixShellName,numStr,'\n')));
    fprintf(fid,char(strcat('#SBATCH --output=',unixShellName,numStr,'.out\n')));
    fprintf(fid,char(strcat('#SBATCH --error=',unixShellName,numStr,'.err\n')));
    fprintf(fid,'#SBATCH --time=2:00:00\n');
    fprintf(fid,'#SBATCH --qos=normal\n');
    fprintf(fid,'#SBATCH --nodes=1\n');
    fprintf(fid,'#SBATCH --mem=10000\n');
    fprintf(fid,'#SBATCH --mail-type=ALL\n');
    fprintf(fid,'#SBATCH --mail-user=ssuen@stanford.edu\n');
    fprintf(fid,'#################\n');
    fprintf(fid,'module load matlab\n');
    fprintf(fid,char(strcat('matlab -nosplash -nodesktop -nodisplay -singleCompThread < $HOME/',sherlockFolderName,'/bar_TBsimulationWrapper',numStr,'.m')));
    fclose(fid);
end


%makes the commands to convert the .scripts to unix
cd(outFoldName)
fileNameStr = 'aaa_convertToUnix.txt';
fid = fopen(fileNameStr,'wt');

for num = totNumScriptsFrom:totNumScriptsTo
    numStr = num2str(num);
    fprintf(fid,char(strcat('dos2unix',{' '}, unixShellName,numStr ,'.script\n')));
    
end
fclose(fid);
