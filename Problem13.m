clc; clear; close all;

% Data Generation
n = 300;

% Piece-wise d_true
n_perm = randperm(n);
piece_number = ceil(n*0.2);
piece_indices = sort(n_perm(1:piece_number));

d_true = zeros(n,1);
prev_ind = 1;
for i=1:piece_number
    next_ind = piece_indices(i);
    d_true(prev_ind:next_ind) = randn(1,1);
    prev_ind = next_ind+1;
end
d_true(prev_ind:n) = rand(1);

% Smooth d_true

% d_true = sign((1:2:n)+rand(1));

d = d_true + randn(n,1)*0.05;
lambda = 0.1;
Dx = @(x) x(1:n-1) - x(2:n);
f = @(x) 0.5*norm(x-d,2)^2 + lambda*norm(Dx(x),1);

%CVX solutions
cvx_begin quiet
    variables x_cvx(n)
    minimize 0.5*pow_pos(norm(x_cvx-d,2),2) + lambda*norm(Dx(x_cvx),1)
cvx_end

% DPG
epsilon = 10^-8;
[x_est_DPG, fk_DPG, iter_DPG] = P13_DPG(d, lambda, f(d_true), epsilon);

% FDPG
[x_est_FDPG, fk_FDPG, iter_FDPG] = P13_FDPG(d, lambda, f(d_true), epsilon);

f0=figure;
stairs(d_true);
hold on
% stairs(x_cvx);
stairs(x_est_DPG);
stairs(x_est_FDPG);
hold off

f1=figure;
semilogy(1:iter_DPG,fk_DPG-f(d_true));
hold on;
semilogy(1:iter_FDPG,fk_FDPG-f(d_true));
hold off;
