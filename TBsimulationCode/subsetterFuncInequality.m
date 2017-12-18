function outIndexVector = subsetterFuncInequality(stateMat, inIndexVector, index, low, high)
% stateMat is the matrix of individual health/status state, 
% inIndexVector is an array of people to consider.  A subset of stateMat, specified by index. 
%        (Use -1 as first element of array if all people should be considered)
% index is the column of the stateMat to be considered, corresponding to sex, age, etc
% low and high are the bounds (inclusive) of permissible values
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

    logicalSubset =   (stateMat(inIndexVector,index) >= low) ...
                    & (stateMat(inIndexVector,index) <= high);

    outIndexVector = inIndexVector(logicalSubset);
end

%%subsetterFunc(stateMat,[-1],%index,%permissible Values)