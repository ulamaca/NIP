function [wSVM, bSVM, a, indexSVs, classSVM] = SVC(xData,xClass,C)
% support vector classification
    % to-be-done1: automatically selecting best C by CV

if nargin<3
    C = 10; % penalization parameters (to tune by CV)
end

N = length(xData);
d = size(xData,2);
k = repmat(xClass,1,d).*xData;
K =  k * k' ;% kernel
f = -ones(N,1);
Aeq = xClass';
beq = 0;
A = [-eye(N); eye(N)];
b = [zeros(N,1); C * ones(N,1)];
% lb = 0; % problem1: using lb <= x <= ub is not working properly...
% ub = C;
%x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
a = quadprog(K,f,A,b,Aeq,beq);

% SVM prediction
wSVM = ( (a.*xClass)' * xData )';
eps = 10^-5;
indexSVs = find(a>eps);% which of input are support vectors
bSVM = xClass(indexSVs(1)) - xData(indexSVs(1),:) * wSVM;
classSVM = 0.5 * sign(xData * wSVM + bSVM) + 1.5 ; % >=1 ; <= -1, what if slack is introduced? (*1 to be correct)

