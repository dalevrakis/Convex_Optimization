function [x_est, f1k_iter, iter] = ISTA(A, b, x_init, f_opt, c, lambda)
    f1 = @(x) (1/2)*norm(A*x-b,2)^2+lambda*norm(x,1);
    ATA = A'*A;
    ATb = A'*b; 
    subg_f = @(x) ATA*x-ATb;
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    iter = 0;
    step = 1/norm(ATA,'fro');
    step_lambda = step*lambda;
    
    %First iteration f_k subg
    f1_k = f1(x_k);
    subg_k = subg_f(x_k);
    
    %
    f1k_iter = zeros(T,1);
    f1k_iter(iter+1) = f1_k;
    for i=1:T-1
        if(f1_k <= c*f_opt)
            break;
        end
        
        x_k = Soft_Thresholding(x_k - step*subg_k, step_lambda);
        
        %Next Iteration f_k subg
        f1_k = f1(x_k);
        subg_k = subg_f(x_k);
        iter = iter + 1;
        
        f1k_iter(iter+1)=f1_k;
    end
    f1k_iter  = f1k_iter(1:iter+1);
    x_est = x_k;
end