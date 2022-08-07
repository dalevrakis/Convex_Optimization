function [x_est, fk_iter, iter] = P9_SFISTA_b(A, b, x_init, f_opt, epsilon)
    F = @(x) norm(x,1);
    P_C = @(x) x - A'*((A*A')\(A*x-b)); %Projection onto Ax-b=0
    n = size(x_init,1);
    
    l_h = sqrt(n);
    mu = epsilon/(l_h^2);
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    y_k = x_init;
    
    iter = 0;
    step = 1;
    
    %First iteration f_k subg
    f_k = F(x_k);
    
    %
    fk_iter = zeros(T,1);
    fk_iter(iter+1) = f_k;
    for i=1:T-1
        if( abs(f_k - f_opt) <= epsilon )
            break;
        end
        
        new_x_k = P_C(Soft_Thresholding(y_k,mu));

        new_step = ( 1+sqrt(1+4*step^2) )/2 ;

        y_k = new_x_k  + ( (step-1)/new_step )*(new_x_k - x_k);
        
        x_k = new_x_k;
        step = new_step;
        
        %Next Iteration f_k subg
        f_k = F(x_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end