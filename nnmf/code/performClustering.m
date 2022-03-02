function clust = performClustering(feat)
% cluster = performClustering(feat,vararing)
%
%   feat:       feature matrix as computed by computeRGCFeatures
%   cluster:    struct containing clustering result
%                   - model
%                   - posterior
%                   - assignments
%                   - bic

disp('FUNCTION: performClustering')
disp('===========================')

p.maxClust = 20;
p.minClust = 1;
p.maxIter = 500;
p.initReplicates = 20;
p.finalReplicates = 20;
p.regularize = 1e-5;
p.covType = 'full';
%p = parseVarArgs(p,varargin{:});

% set random number generator
rng(0);

% fit models with increasing complexity
bic = zeros(p.maxClust-p.minClust+1,1);
opt = statset('MaxIter',p.maxIter);
nClust = p.minClust:p.maxClust;
for j=1:length(nClust)
    gm=gmdistribution.fit(feat,nClust(j),'regularize',p.regularize, ...
        'CovType', p.covType,'Options',opt,'Replicates',p.initReplicates);
    bic(j) = gm.BIC;
    fprintf('Evaluating model with %d clusters. BIC = %.2f \n',nClust(j),bic(j)) 
end

% select best model
[~, mc] = min(bic);

% refit model with optimal parameters
fgm = gmdistribution.fit(feat,nClust(mc),'regularize',p.regularize, ...
    'CovType', p.covType,'Options',opt,'Replicates',p.finalReplicates);

% return structure
clust.model = fgm;
clust.idx = cluster(fgm,feat);
clust.posterior = posterior(fgm,feat);
clust.maxPosterior = max(clust.posterior,[],2);
clust.nClust = nClust;
clust.bic = bic;
clust.K = fgm.NComponents;
clust.nCells = size(feat,1);

disp('')
