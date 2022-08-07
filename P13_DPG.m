function [x_est, fk_iter, iter] = P13_DPG(d, lambda, f_opt, epsilon)
    
    n = size(d,1);
    
    Dx = @(x) x(1:n-1) - x(2:n);
    DTx = @(x) [x ; 0] - [0 ; x];
    f = @(x) 0.5*norm(x-d,2)^2 + lambda*norm(Dx(x),1);
    
    y_k = zeros(n-1,1);
    iter = 0;
    T = 10^8; %Time Horizon
    
    %
    x_k = DTx(y_k) + d;
    y_k = y_k - 0.25*Dx(x_k) + 0.25*Soft_Thresholding(Dx(x_k)-4*y_k,4*lambda);
    
    f_k=f(x_k);
    
    fk_iter = zeros(T,1);
    iter = iter+1;
    fk_iter(iter)=f_k;
    for i=1:T
        if ( f_k-f_opt <= epsilon )
            break;
        end
        
        x_k = DTx(y_k) + d;
        
        y_k = y_k - 0.25*Dx(x_k) + 0.25*Soft_Thresholding(Dx(x_k)-4*y_k,4*lambda);
        
        %Next Iteration f_k
        f_k = f(x_k);
        iter = iter + 1;
        fk_iter(iter)=f_k;
    end
    fk_iter  = fk_iter(1:iter);
    x_est = x_k;
end