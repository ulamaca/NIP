function [] = RasterPlot(tSP)
%% Raster Plot <- spike time can be transformed into more clear ' | ' 
% tSP: N-trial spike trains, in N-d cell array format, each of the entry
% is a vector of time points for a spike train

N = length(tSP);
for i=1:N    
    spikeTrain = tSP{i};
    plot(spikeTrain,i*ones(1,length(spikeTrain)),'*'); hold on
%     for t=1:length(spikeTrain)
%     line([spikeTrain(t) spikeTrain(t)], [i-0.25 i+0.25],'Color','r')
%     end
end
title([num2str(N) ' trials raster plot'])
xlabel('time (sec)')
ylabel('i-th trial')
ylim([0 N+1])    