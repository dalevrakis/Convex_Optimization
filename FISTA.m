function [x_est, f1k_iter, iter] = FISTA(A, b, x_init, f_opt, c, lambda)
    f1 = @(x) (1/2)*norm(A*x-b,2)^2+lambda*norm(x,1);
    ATA = A'*A;
    ATb = A'*b; 
    subg_f = @(x) ATA*x-ATb;
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    y_k = x_init;
    
    iter = 0;
    step = 1;
    L_k = norm(ATA,'fro');
    
    %First iteration f_k subg
    f1_k = f1(x_k);
    subg_k = subg_f(y_k);
    
    %
    f1k_iter = zeros(T,1);
    f1k_iter(iter+1) = f1_k;
    for i=1:T-1
        if(f1_k <= c*f_opt)
            break;
        end
        
        new_x_k = Soft_Thresholding(y_k - (1/L_k)*subg_k, 1/L_k);
        
        new_step = ( 1+sqrt(1+4*step^2) )/2 ;
%         new_step = ( iter+2 )/3 ;
        y_k = new_x_k  + ( (step-1)/new_step )*(new_x_k - x_k);
        
        x_k = new_x_k;
        step = new_step;
        
        %Next Iteration f_k subg
        f1_k = f1(x_k);
        subg_k = subg_f(y_k);
        iter = iter + 1;
        
        f1k_iter(iter+1)=f1_k;
    end
    f1k_iter  = f1k_iter(1:iter+1);
    x_est = x_k;
end