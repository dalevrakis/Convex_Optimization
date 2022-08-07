function [x_est, fk_iter, iter] = P10_Mirror_Descent(A, b, c, x_init, f_opt)
    
    f = @(x) norm(A*x-b,1);
    subg = @(x) A'*sgn(A*x-b);
    m = size(A,1);
    
    x_k = x_init;
    iter = 0;
    T = 10^8; %Time Horizon
    
    %First iteration f_k subg
    f_k = f(x_k);
    subg_k = subg(x_k);
    
    %
    fk_iter = zeros(T,1);
    fk_iter(iter+1) = f_k;
    for i=1:T
        if (f_k <= c*f_opt )
            break;
        end
        
        step = sqrt(2)/( norm(subg_k,Inf)*sqrt(iter+1) );
        
        x_k = x_k.*exp(-step*subg_k)./( sum(x_k.*exp(-step*subg_k)) );
        
        %Next Iteration f_k subg
        f_k = f(x_k);
        subg_k = subg(x_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end