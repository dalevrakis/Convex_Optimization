function [x_est, fk_iter, iter] = P8_SFISTA(A, b, c, x_init, f_opt, lambda, epsilon)
    F = @(x) norm(A*x-b,1)+lambda*norm(x,1);
    [m,n] = size(A);
    ATA = A'*A;
    ATb = A'*b;

    alpha = max(eig(ATA));
    beta = m/2;
    mu = sqrt(alpha/beta)*epsilon/(2*sqrt(alpha*beta));
    subg_Fmu = @(x) (1/mu)*A'*( A*x-b - Soft_Thresholding(A*x-b,mu)) + lambda*sgn(x);
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    y_k = x_init;
    
    iter = 0;
    step = 1;
    L = alpha/mu;
    
    %First iteration f_k subg
    f_k = F(x_k);
    subg_k = subg_Fmu(y_k);
    
    %
    fk_iter = zeros(T,1);
    fk_iter(iter+1) = f_k;
    for i=1:T-1
        if( f_k <= c*f_opt )
            break;
        end
        
        new_x_k = Soft_Thresholding(y_k - (1/L)*subg_k, lambda/L);
%                 new_x_k = y_k - (1/L)*subg_k;

        new_step = ( 1+sqrt(1+4*step^2) )/2 ;

        y_k = new_x_k  + ( (step-1)/new_step )*(new_x_k - x_k);
        
        x_k = new_x_k;
        step = new_step;
        
        %Next Iteration f_k subg
        f_k = F(x_k);
        subg_k = subg_Fmu(y_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end