function [x_est, fk_iter, iter] = P8_Prox_Subgradient(A, b, x_init, f_opt, c ,lambda)

    f = @(x) norm(A*x-b,1) + lambda*norm(x,1);
    subg_f = @(x) A'*sgn(A*x-b);
    m = size(A,1);
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    iter = 0;
    
    %First iteration f_k subg
    f_k = f(x_k);
    subg_k = subg_f(x_k);
    
    %
    fk_iter = zeros(T,1);
    fk_iter(iter+1) = f_k;
    for i=1:T-1
        if(f_k <= c*f_opt)
            break;
        end
        
        step = 1 / (norm(subg_k,2)*sqrt(iter+1));
        
        x_k = Soft_Thresholding(x_k - step*subg_k, lambda*step);
        
        %Next Iteration f_k subg
        f_k = f(x_k);
        subg_k = subg_f(x_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end