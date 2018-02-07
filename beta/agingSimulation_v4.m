%% clear all
close all
clc


%make sampled S, I inputs

%% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exactSol   = 0;  %use approximation to get beta
totSamples = 100;
runStr = 'tb';
if strcmp(runStr, 'gon')
   useGonExample = 1;
else
    useGonExample = 0;
end

% for TB example
inputFold = 'C:\Users\suens\Dropbox\Projects\TB_agePrioritization\simulation\Khyati\june10p5';
inputFold = pwd;
firstAge = 31; %16; %age 15, index 16
lastAge  = 60; %76; %age 75, index 76
ageIntervals = 1;


%for gon SIS example
if useGonExample == 1
    firstAge = 15;
    lastAge = 60;
    inputFold = pwd;
    ageIntervals = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (lastAge - firstAge+1)/ageIntervals;  %total number of age groups.  Start at age 15, go to age 100.

%time horizons
timeHorizon = 25;  %analysis period  ( %how many years into the future to simulate for)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make inputs
% Values needed:
% S, E, mu (natural death), nu (clearance prob), d (TB death prob), N
% (total pop), b (number of new births), alpha_fast, alpha_slow , G
% rows: ages, cols: time.

if useGonExample == 1
    p = makeInputs_v5gon(firstAge, lastAge, inputFold, totSamples,ageIntervals);
else
    p = makeInputs_v5(firstAge, lastAge, inputFold, totSamples,ageIntervals);
end

%make some constants
G = [zeros(1,n);eye(n-1, n)]; %aging matrix
p.T = max(p.dataYrs) - min(p.dataYrs) + 1;

% % fix activation
% multiplier = 0.5; 
% p.alpha_fast = multiplier*p.alpha_fast;
% p.alpha_slow = multiplier*p.alpha_slow;

% SIS model: everyone immediately activates
disp('SIS model: everyone immediately activates, no E compartment')
p.alpha_fast = ones(size(p.alpha_fast));
p.alpha_slow = ones(size(p.alpha_slow));
p.S = p.E + p.S;
p.b = p.b + p.newE;
p.E = zeros(size(p.E));
p.newI = zeros(size(p.newI));
p.newE = zeros(size(p.newE));

%% Test inputs if everyone were healthy (only natural death)
% % multiplier = 1; 0.5; 
% % % [sqS, sqI, sqE, sqN] = ageStratSim(T, G, p.N(:,1), zeros(n,n), zeros(n,1), p.b, zeros(n,T), p.mu, zeros(n,1), zeros(n,1), zeros(n,n), zeros(n,n));
% % % [sqS, sqI, sqE, sqN] = ageStratSim_v2(p.T, G, p.N(:,1),      zeros(n,1) ,zeros(n,1), p.b + p.newE + p.newI, zeros(n,p.T),zeros(n,p.T), p.nu, p.mu, p.d, zeros(n,1), zeros(n,1), zeros(n,n));
% % % [sqS, sqI, sqE, sqN] = ageStratSim_v2(p.T, G, p.S(:,1)+p.I(:,1)+p.E(:,1), zeros(n,1) ,zeros(n,1), p.b + p.newE + p.newI, zeros(n,p.T),zeros(n,p.T), p.nu, p.mu, p.d, zeros(n,1), zeros(n,1), zeros(n,n));
% % % [sqS, sqI, sqE, sqN] = ageStratSim_v2(p.T, G, p.S(:,1)+p.I(:,1), p.E(:,1) ,zeros(n,1), p.b+p.newI,  p.newE, zeros(n,p.T), p.nu, p.mu, p.d, multiplier*p.alpha_fast, multiplier*p.alpha_slow, zeros(n,n));
% [sqS, sqI, sqE, sqN] =   ageStratSim_v2(p.T, G, p.S(:,1), p.E(:,1), p.I(:,1), p.b,  p.newE, p.newI, p.nu, p.mu, p.d, multiplier*p.alpha_fast, multiplier*p.alpha_slow, zeros(n,n) );
% % % [sqS, sqI, sqE, sqN] =   ageStratSim_v2(p.T, G, p.S(:,1), p.E(:,1), p.I(:,1), p.b,  p.newE, p.newI, p.nu, p.mu, p.d, multiplier*p.alpha_fast, multiplier*p.alpha_slow, 0.1*ones(n,n) );
% 
% % set multiplier on activation
% p.alpha_fast = multiplier*p.alpha_fast;
% p.alpha_slow = multiplier*p.alpha_slow;
% 
% % % plot(p.N,sqN); hold on; plot([0 25000000],[0 25000000])
% % plot(sum(p.N),'r'); hold on; plot(sum(sqN),'b');
% makePlots1(p, sqN, sqS, sqE, sqI)
% % % assert(  sum(sum( (abs(p.N-sqN)./p.N) > 0.03)) == 0  )  %no element of the simulated N is off the true N by more than 3%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% THIS WORKS %%%%%%%%%%%%%%%%%%%%
% % %% Calculate feasible beta transmission matrix.
% beta = findFeasibleBeta_v2(exactSol, n, G, p.S, p.I, p.I_dataYr_idx, p.N, p.mu, p.b, p.nu, p.d);
% % %% Run simulation with clearance at status quo levels.  Calibrate activation too
% multiplier = 0.5;
% [simS, simI, simE, simN] = ageStratSim_v2(T, G, p.S(:,1), p.E(:,1) , p.I(:,1), p.b, p.newE, p.newI, p.nu, p.mu, p.d, multiplier*p.alpha_fast, multiplier*p.alpha_slow, beta);
%
%
%
% beta = findFeasibleBeta_v3(exactSol, n, G, p.S, p.I, p.I_dataYr_idx, p.N, p.mu, p.b, p.nu, p.d);
% multiplier = 0.5;
% [simS, simI, simE, simN] = ageStratSim_v2(T, G, p.S(:,1), p.E(:,1) , p.I(:,1), p.b, p.newE, p.newI, p.nu, p.mu, p.d, multiplier*p.alpha_fast, multiplier*p.alpha_slow, beta);
% makePlots1
% %%%%%%%%%%%%%%%%%%  END THIS WORKS %%%%%%%%%%%%%%%%%%

%% Find feasible beta's 
for sampleNum = 1:totSamples
    % Calculate feasible beta transmission matrix from 1993-2006 data.
    beta{sampleNum} = findFeasibleBeta_v4(exactSol, n, G, p.S_sampled{sampleNum}, p.I_sampled{sampleNum}, p.I_dataYr_idx, p.N, p.mu, p.b, p.nu_sampled{sampleNum}, p.d);
    
end   

% % check calibration
% for sampleNum = 1:1
%     %check calibration
%     [cal.S{sampleNum}, cal.I{sampleNum}, cal.E{sampleNum}, cal.N{sampleNum}] ...
%         = ageStratSim_v2(p.T, G, p.S(:,1), p.E(:,1) , p.I(:,1), p.b,  p.newE, p.newI, p.nu, p.mu, p.d, p.alpha_fast, p.alpha_slow, beta{sampleNum} );
% end
% 
% %%% make some plots (should match p)
% makePlots2(p,cal,1)



%% Simulate the SQ future

% Set parameters for simulating into the future
disp('need to set parameters for simulating into the future, multiplier')
% f for future
f.T	=	timeHorizon; %p.T;  %how many years into the future to simulate for
for sampleNum = 1:totSamples
    f.S{sampleNum}	=	[ p.S_sampled{sampleNum}(:,p.T), nan(n,f.T)];  %starting at the last year for which we have data
    f.E{sampleNum}	=	[ p.E(:,p.T), nan(n,f.T)];
    f.I{sampleNum}	=	[ p.I_sampled{sampleNum}(:,p.T), nan(n,f.T)];
    f.nu{sampleNum} =	repmat(p.nu_sampled{sampleNum}(:,p.T),1,f.T);
end

% these stay the same
f.b             =	 repmat(p.b(:,p.T-1),1,f.T);  %use second to last
f.newE      	=	 repmat(p.newE(:,p.T-1),1,f.T);%use second to last
f.newI        	=	 repmat(p.newI(:,p.T-1),1,f.T); %use second to last
f.mu	        =	 repmat(p.mu(:,p.T),1,f.T);
f.d	            =	 p.d(:,1);  %stationary
f.alpha_fast	=	 p.alpha_fast;    %remember, no multipliers!
f.alpha_slow	=	 p.alpha_slow;

% Simluate the future
for sampleNum = 1:totSamples
    % Run simulation with clearance at status quo levels.  Calibrate activation too    
    [f.S{sampleNum}, f.I{sampleNum}, f.E{sampleNum}, f.N{sampleNum}] = ...
        ageStratSim_v2(f.T, G, f.S{sampleNum}(:,1), f.E{sampleNum}(:,1) , f.I{sampleNum}(:,1), f.b, f.newE, f.newI, f.nu{sampleNum}, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta{sampleNum});

end

%make some plots (these shouldn't match, p is just drawn in as reference)
% makePlots2(p,f,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Generate data simulation with an example beta to see if it's recoverable
% % beta = reshape(onePossibleBeta+(nullSpaceofA*constants),n,n);  % DBEUGGING PLACEHOLDER FOR NOW
% beta = 10*eye(n,n);
% for i = 6:n
%     j = i-5;
%     beta(i,j) = 5;
% end
%
% % [simS, simI, simE, simN] = ageStratSim(T, G, p.S(:,1), repmat(p.E(:,1),1,n) , p.I(:,1), p.b, p.nu, p.mu, p.d, p.alpha_fast, p.alpha_slow, beta);
% [simS, simI, simE, simN] = ageStratSim_v2(5, G, p.S(:,1), p.E(:,1) ,p.I(:,1), p.b,  p.newE, p.newI, p.nu, p.mu, p.d, p.alpha_fast, p.alpha_slow, beta);
%
% I_dataYr_idx = 1:size(simS,2)-1;
% % I_dataYr_idx = 1:4;
% testBeta = findFeasibleBeta_v2(exactSol, n, G, simS, simI, I_dataYr_idx, simN, p.mu, p.b, p.nu, p.d);
%
% figure
% imagesc(testBeta); colorbar;  title('Found Beta')
% figure
% imagesc(beta); colorbar;  title('True Beta')
%
% disp('what')



