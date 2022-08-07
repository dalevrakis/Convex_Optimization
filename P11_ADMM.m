function [x_est, fk_iter, iter] = P11_ADMM(A, b, c, rho, x_init, f_opt, e, max_iter)
    
    f = @(x) c'*x;
    [m,n] = size(A);
    AAT = A*A';
    
    x_k = x_init;
    z_k = x_k;
    u_k = x_k;
    
    iter = 0;
    T = 10^8; %Time Horizon
    
    %First iteration f_k subg
    f_k = f(x_k);
    
    %
    fk_iter = zeros(max_iter,1);
    fk_iter(iter+1) = f_k;
    for i=1:max_iter
        if ( abs(f_k-f_opt) <= e )
            break;
        end
        
        diff_zu = z_k-u_k; 
        
        v = AAT\( -A*c + rho*A*diff_zu-rho*b);
        
        x_k =-(1/rho)*c + diff_zu-(1/rho)*A'*v;
        
%         if norm(A*x_k - b,2)>= (10^-5)
%             fprintf("Not in domf\n");
%             disp(norm(A*x_k - b,2))
%         end
        
        z_k = max(x_k+u_k,0);
        
        u_k = u_k + x_k - z_k;
        
        iter = iter + 1;
        
        f_k = f(x_k);
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end