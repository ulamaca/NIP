function [Npc, Um, Y, X_hat, Eigv] = pPCA(X,p) 
% Goal: find 'feature' (linear subspace) that accounts for p of total
% variance for a given data set X
% Inputs: 
%   X: in R^(d by N), d=dimension of data and N=total # of samples
%   p: proporsion of variance to be accounted for
% Outputs:
%   Npc: # of extracted PCs
%   Um: the first Npc Principal components
%   Y:  the dimension reduced data (dim=Npc)
%   X_hat: reconstructed data from Y to original space (dim=d)
%   Eigv:  eigenvalues of covX, which means the variance along the
%   direction of the coresponding eigenvector

N = length(X);
d = size(X,1);
Mu_x = mean(X,2);
Xc = X- Mu_x * ones(1,N);% may uses repmat(Mu_x,N), but it exceed the limit
CovX = (1/N) * Xc*Xc'; % unbiased esti. for covX
[U,D,~] = svd(CovX);
Eigv = diag(D);
Npc = find( (cumsum(Eigv)./sum(Eigv))>p ,1,'first');
Um = U(:,1:Npc);    
Y = Um'*X;    % dimension-reduced data
X_hat = Um*Y;

end