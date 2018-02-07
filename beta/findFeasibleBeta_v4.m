function beta = findFeasibleBeta(exactSol, n, G, S, I, I_dataYr_idx, N, mu, b, nu, d)
disp('making feasible beta matrix')

%elements of sparse I/N matrix
nonZeroElements_x = reshape(repmat([1:n],n,1),[],1);
nonZeroElements_y = [];
for startIdx = 1:n
    nonZeroElements_y = [nonZeroElements_y;(startIdx:n:(n*n))'];
end

%construct equation to solve for beta
yFinal = [];
Afinal = [];

allPeriods = (1:size(S,2)-1);
for t = allPeriods(I_dataYr_idx)
    fi = repmat(I(:,t)./N(:,t),n,1);
    vectorizedIoverN = sparse(nonZeroElements_x,nonZeroElements_y,fi,n,n*n);
    A{t} = -G*diag(S(:,t))*diag(1-mu(:,t))*vectorizedIoverN;
    y{t} = S(:,t+1)-b(:,t)-G*(diag(1-mu(:,t))*S(:,t) + diag(nu(:,t))*diag(1-d)*I(:,t));
    
    %stack time periods
    yFinal = [yFinal; y{t}(2:end,:)];
    Afinal = [Afinal; A{t}(2:end,:)];
end

% % DEBUGGING %does it work at time period 1
% for t = 1
%     onePossibleBeta = pinv(A{t})*y{t};
%     nullSpaceofA = null(A{t});
%
%     %check this worked
%     assert(sum(abs(y{t} - A{t}*onePossibleBeta)) <= 0.000000001);
% end

%solve for beta
% onePossibleBeta = pinv(Afinal)*yFinal;
onePossibleBeta = Afinal\yFinal;

nullSpaceofA = null(Afinal);

%check this worked
assert(sum(abs(yFinal - Afinal*onePossibleBeta)) <= 0.00001);


%get some feasible betas: need to be in the positive orthant   %DEBUGGING NO NEGATIVES
% for iter = 1:totNumBetasSampled
%     constants = repmat(rand(1,size(nullSpaceofA,2)), size(onePossibleBeta,1),1);
%     feasibleBetas(:,iter) = onePossibleBeta + sum(constants.*nullSpaceofA,2);
% end




%% Exact solutions
% rankOfA = size(nullSpaceofA,2);
% % options = optimoptions('linprog','Algorithm','interior-point','Display','iter')
% options = optimoptions('linprog','Algorithm','simplex','Display','iter');
% beta = linprog(ones(rankOfA,1),-nullSpaceofA,-onePossibleBeta,[],[],[],[],[],options);
cvx_solver_settings -clear
cvx_quiet false
cvx_solver_settings('TIMELIMIT', 60)
% cvx_solver_settings('TIMELIMIT', 'None')

%using CVX
if exactSol == 1
rankOfAnull = size(nullSpaceofA,2);
cvx_begin 
    variable constants(rankOfAnull)
    minimize( norm(onePossibleBeta+(nullSpaceofA*constants),2))
    onePossibleBeta+(nullSpaceofA*constants) >= 0
cvx_end

beta = reshape(onePossibleBeta+(nullSpaceofA*constants),n,n);
end

if (exactSol == 0 || isempty(strmatch(cvx_status, 'Solved')))
    disp('Using approximate solution')
    %try the best feasible thing
    cvx_begin 
        variable x(n*n,1)
        minimize( norm(yFinal - (Afinal*x),2) )
        x >= 0;
        x <= 1;
        
        %diagonal is largest value on col
        for i = 1:n 
            for j =((i-1)*n)+1:(i*n); 
                x(((i-1)*n)+i) >= x(j)
            end
        end
        
        %ages are "close" to one another
        for i = 2:n
            for j = 1:n
                colShift = (n*(j-1));
                x(i +colShift) - x(i-1  +colShift) <= 0.03;
                x(i-1 +colShift) - x(i  +colShift) <= 0.03;
            end
        end
        
    cvx_end
    
    beta = reshape(x,n,n);
end

disp('done making feasible beta matrix')