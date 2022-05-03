clc; clear; close all;

% Data Generation
n = 50;
m = 250;

A_m = zeros(2,m,n);
b_m = zeros(m,2);

A_m(1,:,:) = rand(m,n); %Uniform Distribtion
b_m(:,1) = rand(m,1);

A_m(2,:,:) = randn(m,n); %Gaussian Distribusion
b_m(:,2) = randn(m,1);

for dist=1:2
    A = squeeze(A_m(dist,:,:));
    b = b_m(:,dist);
    cvx_begin quiet
        variables x_opt(n)
        minimize norm(A*x_opt-b,1)
    cvx_end

    % if strcmp('Infeasible',cvx_status) == 1
    %     fprintf("Problem Infeasible\n");
    %     return
    % end

%     f = @(x) norm(A*x-b,1);
%     subg = @(x) A'*sgn(A*x-b);
    
    x_0 = rand(n,1);
    c = 1.01;

    %Fixed Polyak Step Size
    [x_est_polyak, fk_polyak, iter_polyak] = Subgradient_Descent(A,b,x_0,cvx_optval,c,0,@f,@subg);
    f_k_best_polyak = f_best_calc(fk_polyak);
    f0 = figure;
    semilogy(1:iter_polyak+1,f_k_best_polyak-cvx_optval);

    %Dynamic Step Size
    [x_est_dyanmic, fk_dynamic, iter_dynamic] = Subgradient_Descent(A,b,x_0,cvx_optval,c,1,@f,@subg); 
    f_k_best_dynamic = f_best_calc(fk_dynamic);
    f1 = figure;
    semilogy(1:iter_dynamic+1,f_k_best_dynamic-cvx_optval);

    %Dynamic-Polyak Comparison
    f2 = figure;
    semilogy(1:iter_polyak+1,f_k_best_polyak-cvx_optval);
    hold on
    semilogy(1:iter_dynamic+1,f_k_best_dynamic-cvx_optval);
    hold off

    %lipschitz constant
    L_f = sqrt(m)*sqrt( max( eig(A.'*A) ) );

    %Upper Bound
    max_iter = max(iter_polyak, iter_dynamic);
    UB = Upper_Bound(max_iter+1,L_f,x_0,x_opt);

    %Dynamic-Polyak-Upper Bound Comparison
    f3 = figure;
    semilogy(1:iter_polyak+1,f_k_best_polyak-cvx_optval);
    hold on
    semilogy(1:iter_dynamic+1,f_k_best_dynamic-cvx_optval);
    semilogy(1:max_iter+1,UB);
    hold off

    %Stochastic Subgradient
    % Theta = 0.5*norm(x_0-x_opt)^2;
    % x_est_stoch = x_0;
    % L_f_stoch = sqrt(m(1))*norm(A,"fro");
    % 
    % f_epoch_stoch = zeros(iter_dynamic+1,1);
    % f_epoch_stoch(1) = f(x_0);
    % for i=1:iter_dynamic
    %    [x_est_stoch, iter_stoch] = Stochastic_Subgradient_Descent(A,b,x_est_stoch,Theta,L_f_stoch);
    %    f_epoch_stoch(i+1) = f(x_est_stoch);
    % end

    [x_est_stoch, f_epoch_stoch] = Stochastic_Subgradient_Descent2(A, b, x_0, x_opt, iter_dynamic);
    f4 = figure;
    semilogy(1:iter_dynamic+1,f_epoch_stoch-cvx_optval);

    %Incremental Subgradient
    [x_est_inc, f_epoch_inc] = Incremental_Subgradient_Descent(A, b, x_0, x_opt, iter_dynamic);
    f5 = figure;
    semilogy(1:iter_dynamic+1,f_epoch_inc-cvx_optval);

    %Complete figures
    f6 = figure;
    semilogy(1:iter_polyak+1,f_k_best_polyak-cvx_optval);
    hold on
    semilogy(1:iter_dynamic+1,f_k_best_dynamic-cvx_optval);
    semilogy(1:max_iter+1,UB);
    semilogy(1:iter_dynamic+1,f_epoch_stoch-cvx_optval);
    semilogy(1:iter_dynamic+1,f_epoch_inc-cvx_optval);
end

function res = f(A,b,x)
    res = norm(A*x-b,1);
end

function res = subg(A,b,x)
    res = A'*sgn(A*x-b);
end