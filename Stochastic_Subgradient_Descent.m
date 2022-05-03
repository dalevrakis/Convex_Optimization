function [x_est, iter] = Stochastic_Subgradient_Descent(A, b, x_init, Theta, L_f)
    
    %Useful variables
    m = size(A,1);
    l = randi(m,m,1);
    
    %Initial Values
    x_k = x_init;
    iter = 0;
    
    for i=1:m
        line = l(i);
        d = sgn( A(line,:)*x_k-b(line) );
%         if d ~= 0
%             subg_ki = d*A(line,:)';
%         else
%             subg_ki = (-1 + 2.*rand(1))*A(line,:)';
%         end
        subg_ki = d*A(line,:)';
        
        step = (sqrt(2*Theta)*m)/(L_f * sqrt(iter+1));
        x_k = x_k - step*subg_ki;
       
        iter = iter + 1;
    end
    x_est = x_k;
end