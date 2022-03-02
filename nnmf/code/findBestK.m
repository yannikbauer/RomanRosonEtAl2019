

%% 1. COPHENETIC CORRELATION COEFFICIENT - PARALLEL

% You repeat NMF several times per rank and you calculate how similar 
% the results are. In other words, how stable the identified clusters are,
% given that the initial seed is random. Choose the highest K before the
% cophenetic coefficient drops.

k_min   = 2;    % fixed value, don't change
k_max   = 50;
maxIter = 200;

option.dis = 0; % sparsennmf info will be suppressed

% assure a random seed
S = RandStream('mt19937ar');
RandStream.setGlobalStream(S);

X = psth.psth;

% preallocate variables
consens_mean = zeros(size(X,1),size(X,1),k_max-k_min+1);
d_res = zeros(k_max-k_min+1,maxIter);

myPool = parpool(4); % 4-cored iMac
parfor ik = k_min:k_max
    
    % pre-allocate variables for parallel processing in nested loops
    consens_iter = zeros(size(X,1),size(X,1),maxIter);
    tmp_d_res = zeros(maxIter,1);
    
    for iter = 1:maxIter % nnmf iteration
        warning off MATLAB:nargchk:deprecated
        
        fprintf('\n... iteration %d; k = %d\n', iter,ik)
        fprintf('=======================\n')
        
        % Sparse-NMF optimized by NNLS
        [A_w,Y_f,~,~,D_r] = sparsenmfnnls(X,ik,option);
        
        % compute connectivity matrix C of size MxM, with entry cij=1 if
        % samples i and j belong to the same cluster, and cij=0 if they belong
        % to different clusters
        [~,clustIdx] = max(A_w,[],2);
        
        conmat = zeros(size(A_w,1));
        pairnums = nchoosek(1:size(conmat,1),2);
        
        for ipair = 1:size(pairnums,1)
            
            pair = pairnums(ipair,:);
            
            var1 = pair(1);
            var2 = pair(2);
            
            % compares the occurrence of the same unit pairs in one feature space
            x = clustIdx(var1) == clustIdx(var2);
            
            conmat(var1,var2) = x;
            conmat(var2,var1) = x;
        end
        consens_iter(:,:,iter) = conmat;
        tmp_d_res(iter) = D_r;
    end
    consens_mean(:,:,ik-1) = mean(consens_iter,3);
    d_res(ik-1,:) = tmp_d_res;
end

delete(myPool)

%% compute cophenetic correlation
cophval = zeros(size(consens_mean,3),1);
for ic = 1:size(consens_mean,3)
    
    % 1. define a consensus matrix for one rank
    C = consens_mean(:,:,ic);
    
    % 2. computer a distance matrix and set 
    C_dist = 1-C;
    C_dist(logical(eye(size(C_dist)))) = 0;
    
    % 3. simulate a pdist output
    P = squareform(C_dist,'tovector');  
    
    % 4. apply linkage
    Z = linkage(P,'average');
    
    % 5. apply cophenet function
    cophval(ic) = cophenet(Z,P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('rev_coph_01h','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot cophenetic correlation
figure; plot(k_min:k_max,cophval,'-o','LineWidth',1.5)
set(gca,'FontSize',11,'TickDir','out','XMinorTick','on','YMinorTick','on')

box off;
xlabel('# of features','FontSize',14);
ylabel('Cophenetic correlation coefficient','FontSize',14)

title('Cophenetic Correlation')



%% 2. RSS against randomized data

% For any dimensionality reduction approach, there is always a loss of
% information compared to your original data (estimated by RSS). Now
% perform NMF for increasing K and calculate RSS with both your original
% dataset and a randomized dataset. When comparing RSS in function of K,
% the RSS decreases with increasing K in the original dataset, but this is
% less the case for the randomized dataset. By comparing both slopes,
% there should be an K where they cross. In other words, how much
% information could you afford to lose (=highest K) before being within
% the noise.

k_min   = 2;    % fixed value, don't change
k_max   = 50;
maxIter = 50;

option.dis = 0; % sparsennmf info will be suppressed

% assure a random seed
S = RandStream('mt19937ar');
RandStream.setGlobalStream(S);

X = psth.psth;       % changed orientation of the matrix

% pre-allocate variables
rsq_ori  = zeros(k_max-k_min+1,maxIter);
rsq_perm = zeros(k_max-k_min+1,maxIter);
d_ori    = zeros(k_max-k_min+1,maxIter);
d_perm   = zeros(k_max-k_min+1,maxIter);

myPool = parpool(4); % 4-cored iMac
parfor ik = k_min:k_max
    
    % pre-allocate variables for parallel processing in nested loops
    tmp_perm   = zeros(maxIter,1);
    tmp_ori    = zeros(maxIter,1);
    tmp_d_perm = zeros(maxIter,1);
    tmp_d_ori  = zeros(maxIter,1);
    
    for iter = 1:maxIter
        warning off MATLAB:nargchk:deprecated
        
        fprintf('\n... iteration %d; k = %d\n', iter,ik)
        fprintf('=======================\n')
        
        permmat = zeros(size(X));
        for icol = 1:size(X,2)
            permmat(:,icol) = randperm(size(X,1));
        end
        
        X_rand = X(permmat);  % randomized dataset
        
        % Original, local minima handling
        [A_ori,Y_ori,~,~,D_res] = sparsenmfnnls(X,ik,option);
               
        R = X-A_ori*Y_ori;
        RSS = sum(R.^2);
        RSS_mean = mean(RSS);
        tmp_ori(iter) = RSS_mean;
        tmp_d_ori(iter) = D_res;
        
        % Permutation
        [A_rand,Y_rand,~,~,D_res] = sparsenmfnnls(X_rand,ik,option);
        
        R = X_rand-A_rand*Y_rand;
        RSS = sum(R.^2);
        RSS_mean = mean(RSS);
        tmp_perm(iter) = RSS_mean;
        tmp_d_perm(iter) = D_res;
    end
    
    rsq_ori(ik-1,:)  = tmp_ori;
    rsq_perm(ik-1,:) = tmp_perm;
    d_ori(ik-1,:)    = tmp_d_ori;
    d_perm(ik-1,:)   = tmp_d_perm;
end

delete(myPool);

% determine best k
rsq_ori_mean = min(rsq_ori,[],2);
rsq_perm_mean = mean(rsq_perm,2);

k_seq = k_min:k_max;
[~,ind] = min(abs(rsq_perm_mean-rsq_ori_mean));
k_best = k_seq(ind);


% As the NMF finds different solutions for different initial conditions,
% the factorizations were repeated 100 times using the previously
% determined rank and evaluated according to their SSE. The A/Y-pair with
% the smallest RE was selected for further analysis.

ntrials = 100;

A_tmp = zeros(size(X,1),k_best,ntrials);
Y_tmp = zeros(k_best,size(X,2),ntrials);

D_res = zeros(ntrials,1);
RSS_valid = zeros(ntrials,1);

myPool = parpool(4); % 4-cored iMac
parfor iter = 1:ntrials    
    warning off MATLAB:nargchk:deprecated
        
    fprintf('\n... iteration %d; k = %d\n', iter,k_best)
    fprintf('=======================\n')
    
    [A,Y,~,~,D] = sparsenmfnnls(X,k_best,option);
    R = X-A*Y;
    RSS = sum(R.^2);
        
    RSS_valid(iter) = mean(RSS);
    A_tmp(:,:,iter) = A;
    Y_tmp(:,:,iter) = Y;
    D_res(iter) = D;
end

delete(myPool);

% select best A/Y
[~,ind] = min(RSS_valid);
A = A_tmp(:,:,ind);
Y = Y_tmp(:,:,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('rev_SSE_01h','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT
figure; hold on
sem_re = std(rsq_ori,[],2)/sqrt(size(rsq_ori,1));
sem_perm = std(rsq_perm,[],2)/sqrt(size(rsq_perm,1));

ftr1 = shadedErrorBar(k_min:k_max,rsq_ori_mean,sem_re,'k','transparent');
set(ftr1.mainLine,'LineWidth',1,'Color','r')

ftr2 = shadedErrorBar(k_min:k_max,rsq_perm_mean,sem_perm,'k','transparent');
set(ftr2.mainLine,'LineWidth',1,'Color','b')

set(gca,'FontSize',11,'TickDir','out','XMinorTick','on','YMinorTick','on')

xlabel('# of features','FontSize',14)
ylabel('Sum of Square Error','FontSize',14)
title('SSE against randomized data')

h = findobj(gca,'Type','line');

lgd = legend(h([1,4]),'data','permutation');
lgd.LineWidth = 0.5;


% plot feature line
yval = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);

xval = zeros(100,1);
xval(1:100,1) = k_best;

plot(xval,yval,'--', 'color', [0.2 0.2 0.2],'LineWidth',1.5)



