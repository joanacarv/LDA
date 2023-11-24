function [v,c] = ldareg(x,y,gamma)

% [v,c] = ldareg(x,y,gamma)
% Implements regularized fisher's linear discriminant analysis
%
% x - N input samples [N,D]
% y - N binary labels {0,1}
% gamma - shrinkage constant [0,1]

[N,D]=size(x);
x=x';

% for each class of trials compute mean and covariance 
tmpi = find(y); N1=length(tmpi);
mX1 = mean(x(:,tmpi),2);
S1 = cov(x(:,tmpi)');
S1 = (1-gamma)*S1 + gamma*(trace(S1)/D)*eye(D);

tmpi = find(y==0); N2=length(tmpi);
mX2 = mean(x(:,tmpi),2);
S2 = cov(x(:,tmpi)');
S2 = (1-gamma)*S2 + gamma*(trace(S2)/D)*eye(D);

% compute common covariance matrices
S = (N1*S1+N2*S2)/N;

lastwarn('');
tmpS = geninv(S);
if strncmp(lastwarn,'Matrix is close to singular or badly scaled.',44)
    tmpS = pinv(S);
end

% compute discriminating vector
v = tmpS*(mX1-mX2);
% compute threshold (to move discriminating criterion at origin)
c = v'*(mX1+mX2)/2;