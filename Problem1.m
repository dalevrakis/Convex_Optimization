clc; clear;
n=2;

e = ones(n,1);
x_0 = randn(n,1)*5;

eq = @(m) e'*max(x_0 - m,0) - 1;

% Bisection

lb = max(x_0)+1; %lower bound
ub = sum(x_0)/n+1; %upper bound
step = 0.1;

MaxIter = 10^5;
err = 10^-8;
mu = 0;
r = 1;

%Find good lower bound
while eq(lb) < 0
    lb = lb - step;
end

%Find good upper bound
while eq(ub) > 0
    ub = ub + step;
end

iter = 0;
for i = 1:MaxIter
    if(r<err)
        break;
    end
    
    %Calc mid
    mu = (ub-lb)/2;
    r = eq(mu);
    if r > 0
        lb = mu;
    else
        ub = mu;
    end
    iter = iter + 1;
end

Projection = max(x_0 - mu,0);
cvx_begin
    variables x(n)
    minimize norm(x-x_0)
    subject to
        e'*x==1;
        x>=0;
cvx_end