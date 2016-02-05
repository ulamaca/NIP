%% 1) Example of applying pPCA 
    % Covariance Estimation and Cumulative Covariance plot
    load('mnist_train.mat')
    X0 = train{10};
    X1 = train{1};
    X = [X0, X1];
    N = length(X); % # of total data
    dim = size(X,1);
    mX = mean(X,2);    
    p = 0.7;
    covX = cov(X');
    [NpcX, Um,  Y, X_hat, D] =pPCA(X,p);    
    cVarX = cumsum(D)./sum(D);
    CovX_hat = cov(X_hat');
    ErrRec = mean( sqrt(sum((X_hat - X).^2,1)) ./ sqrt(sum(X.^2,1))); % reconstruction error, should < p
    %ErrRec1 = mean( norm((X_hat-X),2).^2 ./ norm(X,2).^2)
    
    % Plotting
    figure,
    eps = 0.03; % for texting Npcs
    subplot(3,2,1), imagesc( reshape(mX,sqrt(dim),sqrt(dim)) ), colorbar, axis equal
                    title('Data Mean')
    subplot(3,2,2), imagesc( reshape(diag(covX),sqrt(dim),sqrt(dim)) ), colorbar, axis equal
                    title('Data Variance')
    subplot(3,2,3), imagesc( reshape(mean(X_hat,2),sqrt(dim),sqrt(dim)) ), colorbar, axis equal
                    title('Reconstructed Data Mean')
    subplot(3,2,4), imagesc( reshape(diag(CovX_hat),sqrt(dim),sqrt(dim)) ), colorbar, axis equal
                    title('Reconstructed Data Variance') 
    subplot(3,2,[5 6]), plot(1:dim,cVarX,[NpcX NpcX],[0 cVarX(NpcX)],'*r'), xlabel('n-th dimension')                    
                     line([NpcX NpcX],[0 cVarX(NpcX)],'Color','r'); % first para; range of x, 2nd para: range of y                    
                     text(NpcX,eps,['  m_{' num2str(100*p) '%}=' num2str(NpcX)])
                     title(['#PCs explaining ' num2str(100*p) '% variance, reconstructed error=' num2str(ErrRec)])
    
    
