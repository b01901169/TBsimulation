clear
cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\July9')
%string = textread('landUnits.txt', '%s');
string = textread('folderNames.txt', '%s');
j = 1;
for i = 1:size(string,1)
    if strcmp(string{i},'Directory')
        if strcmp(string{i+1},'of')
            if strcmp(string{i+2}(end-2), '-')
                graphFolder{j} = string(i+2)
                holdThis = graphFolder{j};
                startStr{j} = holdThis{1}(end-23:end-20);
                caseStr{j} = holdThis{1}(53:93);
                discontinuity{j} = str2num(startStr{j});
                j = j+1;
            end
        end
    end
end
for j = 1:size(graphFolder,2)
    graphFolder{j}
end

%dlmwrite('folderNames.txt', outStr);

%fprintf(fileID,formatSpec,A1,...,An)
