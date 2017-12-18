
masterFoldStr = 'C:\Users\ssuen\Dropbox\TBproject\code\outputs\jan14\';
folderList = {'base_60perc\2013-01-11_04-10-10'...
'aveTreat_aveCatIV_empUpta_slowDST_1997\2013-01-11_17-27-31'...
'aveTreat_aveCatIV_empUpta_slowDST_2007\2013-01-11_17-35-33'...
'aveTreat_aveCatIV_empUpta_slowDST_2017\2013-01-11_17-43-05'...
'aveTreat_aveCatIV_empUpta_slowDST_2027\2013-01-11_17-43-05'};


dispPercMDRofNewCases = [];
displayOutMdrMeas = [];


%loop over all folders
for i = 1:size(folderList,2);
folder = strcat(masterFoldStr, folderList{i});

currentDirectory = pwd;
	cd(folder)
		newMDRcases = dlmread('numNewMDR_evolvedOrTrans_postBurnIn.csv', ',' ,1,0);
		incMat = dlmread('incidenceMatrix_forMakingWHOComparison.csv', ',' ,1,0);  
		diagnosedPpl = dlmread('diagnosedPpl.csv', ',' ,1,0); 

	cd(currentDirectory)


    percMDRofNewCases_2008 = sum(sum(newMDRcases(145:156,1:2)))./sum(incMat(1704:1715,1));
    
	incRatNumPerHunThou_2008 = sum(incMat(1704:1715,1)./incMat(1704:1715,3))*100000;
	numMDR_prevUntreated_2008 = sum(diagnosedPpl(1704:1715,9));
	numPpl_prevUntreated_2008 = sum(diagnosedPpl(1704:1715,6));
	numrelapse_2008 = sum(diagnosedPpl(1704:1715,11));
	numMDR_retreatment_2008 = sum(diagnosedPpl(1704:1715,10));
	numPpl_retreatment_2008 = sum(diagnosedPpl(1704:1715,7));
	numMDRrelapsers_2008 = sum(diagnosedPpl(1704:1715,12));
	
	incRatNumPerHunThou_2010 = sum(incMat(1728:1739,1)./incMat(1728:1739,3))*100000;
	numMDR_prevUntreated_2010 = sum(diagnosedPpl(1728:1739,9));
	numPpl_prevUntreated_2010 = sum(diagnosedPpl(1728:1739,6));
	numrelapse_2010 = sum(diagnosedPpl(1728:1739,11));
	numMDR_retreatment_2010 = sum(diagnosedPpl(1728:1739,10));
	numPpl_retreatment_2010 = sum(diagnosedPpl(1728:1739,7));
	numMDRrelapsers_2010 = sum(diagnosedPpl(1728:1739,12));

	dispPercMDRofNewCases = [dispPercMDRofNewCases,percMDRofNewCases_2008];
	disOut_2008 = [incRatNumPerHunThou_2008;numMDR_prevUntreated_2008;numPpl_prevUntreated_2008;numrelapse_2008;numMDR_retreatment_2008;numPpl_retreatment_2008;numMDRrelapsers_2008];
	disOut_2010 = [incRatNumPerHunThou_2010; numMDR_prevUntreated_2010;numPpl_prevUntreated_2010;numrelapse_2010;numMDR_retreatment_2010;numPpl_retreatment_2010;numMDRrelapsers_2010];

	displayOutMdrMeas = [displayOutMdrMeas,[disOut_2008;disOut_2010]];

end

origFormat = get(0, 'format');
format long g

disp('ave frac MDR among new TB cases')
mean(dispPercMDRofNewCases)
disp('standard deviation of frac MDR among new TB cases')
std(dispPercMDRofNewCases)

disp('inputs for MDRmatching_validation.csv in the TBproj folder, for 2008 and 2010 (stacked)')
displayOutMdrMeas

set(0,'format', origFormat);
