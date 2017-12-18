function LEbuilderN = getLEbuilderN(folderNameStr)


% cd('C:\Users\ssuen\Dropbox\TBproject\code\outputs\may6_2014\base\b01_2014-05-08_21-14-10')
dir = pwd;
cd(folderNameStr)
filename = 'TB simulation.txt';
fid = fopen(filename);
tline = fgets(fid);

strCaptureOn = 0;

while ischar(tline)
    if strCaptureOn >= 1
        strCaptureOn = strCaptureOn + 1;
    end
    if strcmp(strtrim(tline),'numPpl =')
        strCaptureOn = 1;
    end
    
    if strCaptureOn == 3
        LEbuilderN = str2double(strtrim(tline));
        strCaptureOn = 0;
    end
    
    tline = fgets(fid);
end

cd(dir)