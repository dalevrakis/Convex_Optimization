function [x_est, f_epoch] = Incremental_Subgradient_Descent(A, b, x_init, x_opt, epochs)
    
    %Useful variables
    f = @(x) norm(A*x-b,1);
    m = size(A,1);

    Theta = 0.5*norm(x_init-x_opt)^2;
    L_f = sqrt(m)*norm(A,"fro");

    f_epoch = zeros(epochs+1,1);
    f_epoch(1) = f(x_init);
    
    %Initial Values
    x_k = x_init;
    for i=1:epochs
        for j=1:m
            d = sgn( A(j,:)*x_k-b(j) );
            subg_kj = d*A(j,:)';

            step = (sqrt(Theta))/(L_f*m*sqrt(i*m+j+1));
            x_k = x_k - step*subg_kj;
        end
        f_epoch(i+1) = f(x_k);
    end
    x_est = x_k;
end