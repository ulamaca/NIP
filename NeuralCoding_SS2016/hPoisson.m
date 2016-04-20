function tSP = hPoisson(tau,tEnd)
% A simple homo-Poisson process (spike train) generator
% input: tau -> spiketime parameter
%        tEnd -> end time
tSP = exprnd(tau);
i=2;
while 1
    isi = exprnd(tau);
    tnew = tSP(i-1) + isi;
    if tnew >= tEnd
        break
    end
    tSP(i) = tnew;
    i = i+1;
end


end