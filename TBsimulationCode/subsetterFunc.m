function outIndexVector = subsetterFunc(stateMat, inIndexVector, index, permisVal)
% stateMat is the matrix of individual health/status state, 
% inIndexVector is an array of people to consider.  A subset of stateMat, specified by index. 
%        (Use -1 as first element of array if all people should be considered)
% index is the column of the stateMat to be considered, corresponding to sex, age, etc
% permisVal is an array of all permissible values for stateMat(inIndexVector, index) as to be included in the output 


%if the subset inIndexArray is empty, return nothing.
if isempty(inIndexVector)
    outIndexVector = [];

%if the subset inIndexArray is not empty
else
    %if user wants all people to be considered inIndexVector is all people
    if (inIndexVector(1) == -1)
        inIndexVector = (1:size(stateMat,1))';
    end

    % assume there is at least one entry in the premisVal    
    logicalSubset = (stateMat(inIndexVector,index) == permisVal(1));

    % Now do it again if necessary
    checkIndex = 2;
    while (checkIndex <= length(permisVal))
        logicalSubset = logicalSubset | (stateMat(inIndexVector,index) == permisVal(checkIndex));
        checkIndex = checkIndex + 1;
    end

    % Convert the logical subset into a list of indices
    outIndexVector = inIndexVector(logicalSubset);
end



%%subsetterFunc(stateMat,[-1],%index,%permissible Values)