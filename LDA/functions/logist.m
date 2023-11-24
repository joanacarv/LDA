function [v,D] = logist(x,y)
% [v,D] = logist(x,y)
% Iterative recurcive least squares algorithm for linear logistic model
%
% x - N input samples [N,D]
% y - N binary labels {0,1}

[N,D]=size(x);
v = zeros(D,1); vth=0;

% combine threshold coputation with weight vector.
x = [x ones(N,1)];
v = [v; vth];

% init termination criteria
vold=ones(size(v)); 
count=0;

lambda = [0.5*eps*ones(1,D) 0]';

% clear warning as we will use it to catch conditioning problems
lastwarn('');

% IRLS for binary classification of experts (bernoulli distr.)
while 1
  vold=v;
  mu = bernoull(1,x*v);   % recompute weights
  w = mu.*(1-mu); 
  e = (y - mu);
  grad = x'*e; % - lambda .* v;
  inc = inv(x'*(repmat(w,1,D+1).*x)+diag(lambda)*eye(D+1)) * grad;

  if strncmp(lastwarn,'Matrix is close to singular or badly scaled.',44)
    warning('Bad conditioning. Suggest to reduce subspace.')
    break;
  end

  % update
  v = v + inc; 
 
  % exit if converged
  if norm(vold) & subspace(v,vold)<10^-10, break, end;

  % exit if its taking to long 
  count=count+1;
  if count>100, 
    warning('Not converged after 100 iterations.'); 
    %plot(vnorm)
    break; 
  end;   

  if count==1,  vnorm(count) = NaN; 
  else          vnorm(count) = subspace(v,vold); end;
  
end; % while loop

function [p]=bernoull(x,eta);
% [p] = bernoull(x,eta)
%
% Computes Bernoulli distribution of x for "natural parameter" eta.
% The mean m of a Bernoulli distributions relates to eta as,
% m = exp(eta)/(1+exp(eta));

p = exp(eta.*x - log(1+exp(eta)));








