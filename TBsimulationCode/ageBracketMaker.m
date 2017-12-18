%ageBracketMaker function
function ageBrac = ageBracketMaker(ageRangMat, age)

%initialize lookup table
ageBracLookup = -1*ones(110,1,'int16');

%make lookup table
ageRangMat(ageRangMat == Inf) = 109;
for ii = 1:size(ageRangMat,1)
    ageBracLookup(ageRangMat(ii,1)+1:ageRangMat(ii,2)+1) = ii;
end

%perform lookup
ageBrac = ageBracLookup(age+1);
