clc; clear; close all;

% Data Generation
m = 300;
n = 100;

A = rand(m,n);
b = rand(m,1);

%CVX solutions
lambda = 1;
cvx_begin quiet
    variables x_cvx(n) 
    minimize norm(A*x_cvx-b,1) + lambda*norm(x_cvx,1)
cvx_end

x_0 = 0;
c = 1.1;
%Fixed Polyak Step Size
[x_est_polyak, fk_polyak, iter_polyak] = P8_Subgradient_Descent(A,b,x_0,cvx_optval,c,0,lambda);

%Dynamic Step Size
[x_est_dynamic, fk_dynamic, iter_dynamic] = P8_Subgradient_Descent(A,b,x_0,cvx_optval,c,1,lambda);

%Proximal Subgradient
[x_est_prox, fk_prox, iter_prox] = P8_Prox_Subgradient(A, b, x_0, cvx_optval, c ,lambda);

%SFISTA
epsilon = cvx_optval*(c-1);
[x_est_SFISTA, fk_SFISTA, iter_SFISTA] = P8_SFISTA(A, b, c, x_0, cvx_optval, lambda, epsilon);

%SFISTA smooth g
[x_est_SFISTA_d, fk_SFISTA_d, iter_SFISTA_d] = P8_SFISTA_d(A, b, c, x_0, cvx_optval, lambda, epsilon);

f0=figure;
stairs(x_cvx);
hold on
stairs(x_est_polyak);
stairs(x_est_dynamic);
stairs(x_est_prox);
stairs(x_est_SFISTA);
stairs(x_est_SFISTA_d);
hold off

f1=figure;
semilogy(1:iter_polyak+1,fk_polyak-cvx_optval);
hold on
semilogy(1:iter_dynamic+1,fk_dynamic-cvx_optval);
semilogy(1:iter_prox+1,fk_prox-cvx_optval);
semilogy(1:iter_SFISTA+1,fk_SFISTA-cvx_optval);
semilogy(1:iter_SFISTA_d+1,fk_SFISTA_d-cvx_optval);
hold off
