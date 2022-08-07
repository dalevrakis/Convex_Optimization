clc; clear; close all;

% Data Generation
m = 100;
n = 200;

A = rand(m,n);
x_s = rand(n,1);
b = A*x_s;

%CVX solutions
cvx_begin quiet
    variables x_cvx(n) 
    minimize norm(x_cvx,1)
    subject to
        A*x_cvx-b==0;
cvx_end

x_0 = zeros(n,1);

%SFISTA
epsilon = 10^-2;
[x_est_SFISTA, fk_SFISTA, iter_SFISTA] = P9_SFISTA_b(A, b, x_0, cvx_optval, epsilon);

%PSD polyak step
[x_est_PSD_polyak, fk_PSD_polyak, iter_PSD_polyak] = P9_Subgradient_Descent(A, b, x_0, cvx_optval, 0, iter_SFISTA);

%PSD dynamic step
[x_est_PSD_dynamic, fk_PSD_dynamic, iter_PSD_dynamic] = P9_Subgradient_Descent(A, b, x_0, cvx_optval, 1, iter_SFISTA);

f0=figure;
stairs(x_cvx);
hold on
stairs(x_est_SFISTA);
stairs(x_est_PSD_polyak);
stairs(x_est_PSD_dynamic);
hold off

f1=figure;
semilogy(1:iter_SFISTA+1,fk_SFISTA-cvx_optval);
hold on
semilogy(1:iter_PSD_polyak+1,fk_PSD_polyak-cvx_optval);
semilogy(1:iter_PSD_dynamic+1,fk_PSD_dynamic-cvx_optval);
hold off