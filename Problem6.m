clc; clear; close all;

% Data Generation
n = 450;
s = 10;
m = ceil(2*s*log(n));

%Non zero elements
n_perm = randperm(n);
s_ind = n_perm(1:s);
non_zero_vals = randn(s,1);

%Sparse vector
x_s = sparse(s_ind,1,non_zero_vals,n,1);

%A,b generation
A = rand(m,n);
b = A*x_s;

%CVX solutions
lambda = 1;
cvx_begin quiet
    variables x1_opt(n)
    minimize ((1/2)*pow_pos(norm(A*x1_opt-b,2),2) + lambda*norm(x1_opt,1))
cvx_end
% norm(x_s - x1_opt)

cvx_begin quiet
    variables x2_opt(n)
    minimize ((1/2)*pow_pos( norm(A*x2_opt-b,2),2 ) + lambda*pow_pos( norm(x2_opt,2),2 ) )
cvx_end
% norm(x_s - x2_opt)

f0 = figure;
stairs(x_s)
hold on
stairs(x1_opt)
stairs(x2_opt)
hold off
legend({'$x_s$','$x_1^{opt}$','$x_2^{opt}$'},'Interpreter','latex')
hold on