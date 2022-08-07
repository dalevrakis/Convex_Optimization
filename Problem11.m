clc; clear; close all;

%Data Generation
m = 50;
n = 500;

A = rand(m,n);
x_s = rand(n,1);
b = A*x_s;
c = randn(n,1)*5;
c_T = c';

% cvx_begin quiet
%     variables x_0(n)
%     minimize 0
%     subject to
%         A*x_0==b;
% cvx_end

%CVX solutions
cvx_begin quiet
    variables x_cvx(n)
    minimize c_T*x_cvx
    subject to
        A*x_cvx==b;
        x_cvx>=0;
cvx_end

% ADMM
rho = 10;
max_iter = 10^7;
e = 10^-5;
x_0 = A\b;
[x_est_ADMM, fk_ADMM, iter_ADMM] = P11_ADMM(A, b, c, rho, x_0, cvx_optval, e, max_iter);

f0=figure;
stairs(x_cvx);
hold on
stairs(x_est_ADMM);
hold off

f1=figure;
semilogy(1:iter_ADMM+1,abs(fk_ADMM-cvx_optval));
