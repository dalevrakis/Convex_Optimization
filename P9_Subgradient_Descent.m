function [x_est, fk_iter, iter] = P9_Subgradient_Descent(A, b, x_init, f_opt, step_size_type, max_iter)

    f = @(x) norm(x,1);
    subg = @(x) sgn(x);
    P_C = @(x) x - A'*((A*A')\(A*x-b)); %Projection onto Ax-b=0
    m = size(A,1);
    

    x_k = x_init;
    iter = 0;
    
    %First iteration f_k subg
    f_k = f(x_k);
    subg_k = subg(x_k);
    
    %
    fk_iter = zeros(max_iter,1);
    fk_iter(iter+1) = f_k;
    for i=1:max_iter
        
        if step_size_type == 0
            step = (f_k-f_opt)/( norm(subg_k(:,1))^2 );
        elseif step_size_type == 1
            step = 1/( norm(subg_k(:,1))*sqrt(iter+1) );
        else
            fprintf("Invalid Option for step_size_type\n");
        end
        
        x_k = P_C(x_k - step*subg_k(:,1));
        
        %Next Iteration f_k subg
        f_k = f(x_k);
        subg_k = subg(x_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end