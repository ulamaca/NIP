function [mu, Sigma, priors, df, assignments, loglike, Npara] = sortSpikes(b,K)
% Do spike sorting by fitting a Gaussian or Student's t mixture model.
%   [mu, Sigma, priors, df, assignments] = sortSpikes(b) fits a mixture
%   model using the features in b, which is a matrix of size #spikes x
%   #features. The df output indicates the type of model (Gaussian: df =
%   inf, t: df < inf)
%   Input K
%
%   The outputs are (K = #components, D = #dimensions):
%       mu              cluster means (K-by-D), each 1-by-D
%       Sigma           cluster covariances (D-by-D-by-K), each D-by-D
%       priors          cluster priors (1-by-K)
%       df              degrees of freedom in the model (Gaussian: inf)
%       assignments     cluster index for each spike (1 to K)
%       K               # of cluster
%       loglike         loglikelihood of the final fitted distribution
%       Npara           # of parameters used in the model

df = inf;   % you don't need to use this variable unless you want to 
            % implement a mixture of t-distributions
            
Niter = 1000; % # of iterations
N = length(b);
d = size(b,2);

%% Initialization
posLatent = zeros(K,N); % K-by-N responsibility of each data point (i.e. posterior of latent var)
priors    = ones(1,K)/K;
assignments = zeros(1,N);
Sigma = zeros(d,d,K);
    for k=1:K
    Sigma(:,:,k) = cov(b)./K; % Sigma initialized by sample var divided by K
                              % note1: input should be in d-by-N format!
    end
c = 1.5; % coefficeint for mean initialization                   
mu    = ones(K,1) * mean(b) + rand(K,d) .* (ones(K,1) * c * median(b)); % initia mean by some random fluactuation;

Npara = K * (1 + d + d*(d-1)/2); % # of parameter used
loglike = 0; % log-likelihood
sloglike = zeros(K,N); % store for later loglike computation
%% EM-Algorithms
for iter = 1:Niter
    % E-step    
    for k=1:K
        likelihood = mvnpdf(b,mu(k,:),Sigma(:,:,k))'; % p(xn|zn), 1-by-N
        joint      = likelihood .* priors(k);          % p(xn,zn), 1-by-N
        postLatent(k,:) = joint;                      % p(zn|xn, parameters...)        
        if iter == Niter
            sloglike(k,:) = likelihood;
        end
    end
    postLatent = postLatent ./ (ones(k,1) * sum(postLatent) );% normalizing postLatent
    
    % M-step
    Nk = sum(postLatent,2);% (k-by-1) equivalent samples of k-th cluster
    priors = Nk'./N;        % (1-by-k)
    mu     = (postLatent * b) ./ (Nk * ones(1,d));
    for k=1:K
        bc = b - ones(N,1) * mu(k,:); % centralized data-b
        temp = zeros(d);
        for n=1:N
        temp = temp + postLatent(k,n) * bc(n,:)' * bc(n,:) ;        
        end
        Sigma(:,:,k) = temp/Nk(k);
    end
    
    if iter == Niter
        [~, assignments] = max(postLatent);
        loglike = sum(log(priors * sloglike));
    end
    % log-Likelihood computation
    
end
end