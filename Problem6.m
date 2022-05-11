clc; clear; close all;

% Data Generation
n = 1000;
s = 50;
m = ceil(2*s*log(n));

%Non zero elements
n_perm = randperm(n);
s_ind = n_perm(1:s);
non_zero_vals = randn(s,1);

%Sparse vector
x_s = sparse(s_ind,1,non_zero_vals,n,1);

%A,b generation
A = randn(m,n);
b = A*x_s;

%CVX solutions
lambda = 1;
cvx_begin quiet
    variables x1_opt(n) 
    minimize ((1/2)*pow_pos(norm(A*x1_opt-b,2),2) + lambda*norm(x1_opt,1))
cvx_end
% norm(x_s - x1_opt)
% f1_optval = cvx_optval;
f1 = @(x) (1/2)*norm(A*x-b,2)^2+lambda*norm(x,1);
f1_optval = f1(x_s);

cvx_begin quiet
    variables x2_opt(n)
    minimize ((1/2)*pow_pos( norm(A*x2_opt-b,2),2 ) + lambda*pow_pos( norm(x2_opt,2),2 ) )
cvx_end
% norm(x_s - x2_opt)
f2_optval = cvx_optval;

f0 = figure;
stairs(x_s)
hold on
stairs(x1_opt)
stairs(x2_opt)
hold off
legend({'$x_s$','$x_1^{opt}$','$x_2^{opt}$'},'Interpreter','latex');

c = 1.01;
x_init = zeros(n,1);
[x_est_ISTA, f1k_iter_ISTA, iter_ISTA] = ISTA(A, b, x_init, f1_optval, c, lambda);

f1 = figure;
semilogy(1:iter_ISTA+1,f1k_iter_ISTA-f1_optval);
hold on

[x_est_FISTA, f1k_iter_FISTA, iter_FISTA] = FISTA(A, b, x_init, f1_optval, c, lambda);
semilogy(1:iter_FISTA+1,f1k_iter_FISTA-f1_optval);

[x_est_VFISTA, f1k_iter_VFISTA, iter_VFISTA] = VFISTA(A, b, x_init, f1_optval, c, lambda);
semilogy(1:iter_VFISTA+1,f1k_iter_VFISTA-f1_optval);
legend({'$ISTA$','$FISTA$','$VFISTA$'},'Interpreter','latex');
hold off

f2 = figure;
stairs(x_s)
hold on
stairs(x_est_ISTA)
stairs(x_est_FISTA)
stairs(x_est_VFISTA)
hold off
legend({'$x_s$','$x_{est} ISTA$','$x_{est} FISTA$','$x_{est} VFISTA$'},'Interpreter','latex');
