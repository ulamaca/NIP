function H = entropy_comp(p,q)
% compute entropy dealing with 0 * log(0) problem
% Inputs:
%   p: is the probability mass to compute the entropy (col. vector)
%   q: density q, in case that this density is different from p
% Output:
%   H: the entropy result

% Note this code allows matrix input to compute multiple entropy
%      simultaneously, the form of the matrix p should be
%      p = [p1, p2, p3, ...pn] each pi is a probability mas with K possible values
%             K-by-n matrix
%      

if isrow(p) % if input is a row-vector
    p = p';
end

pM = p;
pM(pM==0) = pM(pM==0) + 1/2; % s.t. log(0) term ~= -Inf

if nargin<2
    H = - sum(p  .* log2(pM), 1)';
elseif nargin==2
    H = - sum(q .* log2(pM), 1)';
end
end