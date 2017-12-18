function le = calculateLifeExpectancy(ages, deaths)
% le = calculateLifeExpectancy(ages, deaths)
% Calculates life expectancy given a vector of ages and a vector of # of
% people who died at that age.
% 
% Make sure ages and deaths are the same dimensions (and that they are
% either both column vectors or both row vectors)
deathFrac = deaths / sum(deaths);
le = sum(ages .* deathFrac);