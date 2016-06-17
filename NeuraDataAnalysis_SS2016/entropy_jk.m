function H = entropy_jk(sampHist)
% function H = entropy_jk(p)
%   sampHist   count for each state
%   H   jackknife estimate of entropy

sampHist = sampHist(sampHist>0);
K = length(sampHist); % dimension of state space
N = sum(sampHist); % data size
pML = sampHist/N; %p_h\i

% Approach1: For-loop
%   d=1/(N-1);
%   pMLi = sampHist/(N-1);
%   Hh_i = zeros(K,1);
%     for k=1:K
%         temp = pMLi; % temp prob for computing Hh_i
%         temp(k) = temp(k) - d;
%         Hh_i(k) = entropy_comp(temp);
%     end
    

% Approach2: Matrix
ph_k = 1/(N-1) * ( repmat(sampHist,1,K) - eye(K) ); % a matrix
Hh_i = entropy_comp(ph_k);        

% Result Computation
mHh_i = Hh_i' * pML;
H = N * entropy_comp(pML) - (N-1) * mHh_i;

% jack-knifed entropy estimator (paninski, p. 1198)



