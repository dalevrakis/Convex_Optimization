clc; clear; close all;

% Data Generation
m = 10;
n = 100;

A = randn(m,n);
b = randn(m,1);
e = ones(n,1);

%CVX solutions
cvx_begin quiet
    variables x_cvx(n)
    minimize norm(A*x_cvx-b,1)
    subject to
        e'*x_cvx==1;
        x_cvx>=0;
cvx_end
% fprintf("CVX ended\n");

x_0 = e/n;
c = 1.01;

% Polyak PSD
[x_est_PSD_polyak, fk_PSD_polyak, iter_PSD_polyak] = P10_Subgradient_Descent(A, b, c, x_0, cvx_optval, 0);

% Dynamic PSD
[x_est_PSD_dynamic, fk_PSD_dynamic, iter_PSD_dynamic] = P10_Subgradient_Descent(A, b, c, x_0, cvx_optval, 1);

% Mirror Descend
[x_est_MD, fk_MD, iter_MD] = P10_Mirror_Descent(A, b, c, x_0, cvx_optval);

f0=figure;
stairs(x_cvx);
hold on
stairs(x_est_PSD_polyak);
stairs(x_est_PSD_dynamic);
stairs(x_est_MD);
hold off

f1=figure;
semilogy(1:iter_PSD_polyak+1,fk_PSD_polyak-cvx_optval);
hold on
semilogy(1:iter_PSD_dynamic+1,fk_PSD_dynamic-cvx_optval);
semilogy(1:iter_MD+1,fk_MD-cvx_optval);
hold off