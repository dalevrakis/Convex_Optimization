clc; clear; close all;

% Data Generation
n = 300;
m = 50;

n_perm = randperm(n);
piece_number = floor(n/20);
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

e = rand(m,1)*0.001;

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
    minimize ((1/2)*pow_pos( norm(A*x_rls-b,2),2 ) + norm(D*x_rls,1))
cvx_end

F = @(x) (1/2)*norm(A*x-b,2)^2+ norm(D*x,1);
f_opt = F(x_pwc);
x_init = rand(n,1);
epsilon = 10^-4;
[x_est_SFISTA, fk_iter_SFISTA, iter_SFISTA] = SFISTA(A, b, D, x_init, f_opt, epsilon);

f0=figure;
semilogy(1:iter_SFISTA+1,fk_iter_SFISTA-f_opt);

f1=figure;
stairs(x_pwc);
hold on
stairs(x_est_SFISTA);
hold off