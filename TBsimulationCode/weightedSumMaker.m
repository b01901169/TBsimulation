function outcome = weightedSumMaker(prevalence, default, ageBracTable, age, sex)

% Initialize everybody to false
outcome = 0;

%initialize lookup table
prevLookup = default*ones(110,size(prevalence,2));
%make lookup table
ageBracTable(ageBracTable == Inf) = 109;
for ii = 1:size(ageBracTable,1)
    prevLookup(ageBracTable(ii,1)+1:ageBracTable(ii,2)+1,1) = prevalence(ii,1);
    if size(prevalence,2) > 1
        prevLookup(ageBracTable(ii,1)+1:ageBracTable(ii,2)+1,2) = prevalence(ii,2);
    end
end

%stack lookup table for linear index
prevLookup = [prevLookup;prevLookup]; 
if size(prevalence,2) > 1
    %stack lookup table for linear index
    prevLookup = [prevLookup(:,1);prevLookup(:,2)];
end

% Find counts
%get linear index for each person
pplIndex = (age+1) + (sex-1)*110;

%find how many ppl in that category
outcome = sum(prevLookup(pplIndex));
% Done!
