clc; clear; close all;

% Data Generation
m = 1000;
n = 300;

n_perm = randperm(n);
piece_number = ceil(n*0.2);
piece_indices = sort(n_perm(1:piece_number));

x_pwc = zeros(n,1);
prev_ind = 1;
for i=1:piece_number
    next_ind = piece_indices(i);
    x_pwc(prev_ind:next_ind) = rand(1);
    prev_ind = next_ind+1;
end
x_pwc(prev_ind:n) = rand(1);
% stairs(x_pwc);

A = rand(m,n);
% A = randn(m,n);

e = randn(m,1)*0.1;

b = A*x_pwc + e;

cvx_begin quiet
    variables x_sls(n) 
    minimize ((1/2)*pow_pos( norm(A*x_sls-b,2),2 ))
cvx_end

ro = 1;
ones_matrix = ones(n,1);
nones_matrix = -ones(n-1,1);
D = ro*(diag(ones_matrix)+diag(nones_matrix,1));

cvx_begin quiet
    variables x_rls(n) 
    minimize (0.5*sum_square_abs(A*x_rls-b) + norm(D*x_rls,1));
cvx_end

F = @(x) 0.5*sum_square_abs(A*x-b) + norm(D*x,1);
f_opt = F(x_pwc);
x_init = zeros(n,1);
epsilon = 10^-3;
[x_est_SFISTA, fk_iter_SFISTA, iter_SFISTA] = SFISTA(A, b, ro, D, x_init, f_opt, epsilon);

[x_est_subg_polyak, fk_iter_subg_polyak, iter_subg_polyak] = RLS_Subgradient_Descent(A, b, D, x_init, f_opt, 0, epsilon);

[x_est_subg_dynamic, fk_iter_subg_dynamic, iter_subg_dynamic] = RLS_Subgradient_Descent(A, b, D, x_init, f_opt, 1, epsilon);

f0=figure;
semilogy(1:iter_SFISTA+1,fk_iter_SFISTA-f_opt);
hold on
semilogy(1:iter_subg_polyak+1,fk_iter_subg_polyak-f_opt);
semilogy(1:iter_subg_dynamic+1,fk_iter_subg_dynamic-f_opt);
hold off

f1=figure;
stairs(x_pwc);
hold on
stairs(x_est_SFISTA);
stairs(x_est_subg_polyak);
stairs(x_est_subg_dynamic);
stairs(x_rls);
hold off
legend({'$x_{pwc}$','$x^{RLS}_{FISTA}$', '$x^{RLS}_{SUBG}-polyak$', '$x^{RLS}_{SUBG}-dynamic$','$x^{RLS}_{CVX}$'},'Interpreter','latex');