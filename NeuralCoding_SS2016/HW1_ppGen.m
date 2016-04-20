% N Point process generation and Analysis
tau = [1 0.1 0.01];

%% Poisson Spike trainsGeneration with Spike Count, Mean and Variance
Np = length(tau); % # of parameters
N = 100; % # of spike trains generation
tSP = cell(Np,N);
nSP = zeros(Np,N);
for i=1:Np
    for j=1:N
    temp = hPoisson(tau(i),10);
    tSP{i,j} = temp;
    nSP(i,j) = length(temp);
    end
end
meanNSP = mean(nSP,2);
varNSP  = var(nSP');

%% Bernoulli Spikes
lam = [0.1 0.01 0.001];
Nl = length(lam);
T = 10; % duration
tSPb = cell(Nl,T);
nSPb = zeros(Nl,T);
for i=1:Nl
    for j=1:N
    temp = rand(1,T) < (1 - lam(i)); % Bernouli Spike Train
    tSPb{i,j} = find(temp==1);
    nSPb(i,j) = sum(temp);
    end
end
meanNSPb = mean(nSPb,2);
varNSPb = var(nSPb');
%%  Raster Plot
h1=figure(1),
spikeTrainPlot = cell(1,3);
for i=1:Np
    K = randi(N);
    spikeTrainPlot(i) = tSP(i,K);    
end
subplot(1,2,1),RasterPlot(spikeTrainPlot) % lam = [1 0.1 0.01] <- the tau corresponding to each raster
legend(['\tau =' num2str(tau(1))], ['\tau =' num2str(tau(2))], ['\tau =' num2str(tau(3))])
title('Poisson Spike Train with 3 parameters' )

spikeTrainPlot = cell(1,3);
for i=1:Np
    K = randi(N);
    spikeTrainPlot(i) = tSPb(i,K);    
end
subplot(1,2,2), RasterPlot(spikeTrainPlot) % lam = [0.1 0.01 0.001] <- the lambda corresponding to each raster
legend(['\lambda =' num2str(lam(1))], ['\lambda =' num2str(lam(2))], ['\lambda =' num2str(lam(3))])
title('Bernoulli Spike Train with 3 parameters' )
saveas(h1,'fig1_raster_plots','bmp')
%% Poisson Histogram
% 1. what should it look like?
% 2. how to put mean and variance together with the hist?
% 3. how to change the size of bin?
h2=figure(2),
for p=1:size(nSP,1)
    subplot(1,3,p),
    [counts,centers] = hist(nSP(p,:)); 
    rat = counts./sum(counts);
    centersRdd = round(centers);
    bar(centers,rat); hold on
    lambda = 1/tau(p);
    % () Poisson fitting
    fitP = poisspdf(centersRdd,10*lambda);
    fitP = fitP./sum(fitP);
    plot(centersRdd,fitP,'r');
    plot(meanNSP(p),0,'g*')
    plot([meanNSP(p)-sqrt(varNSP(p)), meanNSP(p)+sqrt(varNSP(p))], [0 0],'--*y')
    Fano = varNSP(p)/meanNSP(p);
    legend('spike-count ratio','Poisson fitting')
    title(['Spike Count histogram (Poisson), \tau = ' num2str(tau(p))])
    xlabel('# of spikes')
end
saveas(h2,'fig2_histogram_Poisson','bmp')
%% Bernouli Histogram
h3=figure(3),
for p=1:size(nSPb,1)
    subplot(1,3,p),
    [counts,centers] = hist(nSPb(p,:));
    rat = counts./sum(counts);
    bar(centers,rat);  hold on       
    plot(meanNSPb(p),0,'g*')
    plot([meanNSPb(p)-sqrt(varNSPb(p)), meanNSPb(p)+sqrt(varNSPb(p))], [0 0],'--*y')
    Fano = varNSPb(p)/meanNSPb(p);
    legend('spike-count ratio')
    title(['Spike Count histogram (Bernouli), \lambda = ' num2str(lam(p))])
    xlabel('# of spikes')
end
saveas(h3,'fig3_histogram_Bernoulli','bmp')