clc; clear; close all;
n=2;

e = ones(n,1);
x_0 = randn(n,1)*5;

f = @(m) e'*max(x_0 - m,0) - 1;

%CVX Solution
cvx_begin quiet
    variables x_cvx(n)
    minimize 0.5*pow_pos(norm(x_cvx-x_0,2),2)
    subject to
        e'*x_cvx==1;
        x_cvx>=0;
cvx_end

% Bisection

ub = max(x_0)+1; %upper bound
lb = (sum(x_0)-1)/n - 1; %lower bound
step = 0.1;

MaxIter = 10^6;
TOL = 10^-5;
err = 10^-4;
mu = 0;
f_mu = 1;

% % Find good lower bound
% while eq(lb) < 0
%     lb = lb - step;
% end
% 
% % Find good upper bound
% while eq(ub) > 0
%     ub = ub + step;
% end

iter = 0;
x_err = zeros(MaxIter,1);
for i = 1:MaxIter
    if(abs(f_mu)<=err || (ub-lb)/2<=TOL)
        break;
    end
    
    %Calc mid
    mu = (ub+lb)/2;
    f_mu = f(mu);
    if f_mu > 0
        lb = mu;
    else
        ub = mu;
    end
    iter = iter + 1;
    Projection = max(x_0 - mu,0);
    x_err(iter,1) = mse(Projection,x_cvx);
end
x_err = [mse(x_0,x_cvx) ; x_err(1:iter,1)];

Projection = max(x_0 - mu,0);

f0 = figure();
semilogy(0:iter,x_err);

ylabel('\fontsize{18} MSE(x_{bis},x_{cvx})','interpreter','tex');
xlabel('\fontsize{18} Iterations','interpreter','tex');

saveas(f0,fullfile('D:\Documents\Tuc\HMMY\10th Semester\ConvexOptimization','P1-fig0.png'));