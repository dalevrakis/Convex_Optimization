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
        minimize max(abs(A*x_opt-b))
    cvx_end
    
    x_0 = rand(n,1);
    c = 1.01;

    %Fixed Polyak Step Size
    [x_est_polyak, fk_polyak, iter_polyak] = Subgradient_Descent(A,b,x_0,cvx_optval,c,0,@f,@L_Inf_Norm_Subdiff);
    f_k_best_polyak = f_best_calc(fk_polyak);
    f0 = figure;
    semilogy(1:iter_polyak+1,f_k_best_polyak-cvx_optval);
    
    %Dynamic Step Size
    [x_est_dyanmic, fk_dynamic, iter_dynamic] = Subgradient_Descent(A,b,x_0,cvx_optval,c,1,@f,@L_Inf_Norm_Subdiff); 
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
    L_f=0;
    for i=1:m
        L_f = L_f+norm(A(i,:),2);
    end

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
end

function res = f(A,b,x)
    res = max(abs(A*x-b));
end
