function [x_est, fk_iter, iter] = SFISTA(A, b, ro, D, x_init, f_opt, epsilon)
    F = @(x) 0.5*sum_square_abs(A*x-b) + norm(D*x,1);
    n = size(A,2);
    ATA = A'*A;
    ATb = A'*b;
%     ATA_fro = norm(A'*A,'fro');
%     A_fro = norm(A,'fro');
    
%     D_fro = norm(D,'fro');
    DTD = D'*D;
    
    alpha = max(eig(DTD));
    beta = n/2;
    L_f = max(eig(ATA));
    mu = sqrt(alpha/beta)*epsilon/( sqrt(alpha*beta)+sqrt(alpha*beta + L_f*epsilon) );
%     mu = ( (2*D_fro)/sqrt(m) )*epsilon/( sqrt(D_fro^2*m) + sqrt(D_fro^2*m + 2*ATA_fro*epsilon)); 
    
    Dx = @(x)ro*[x(1:end-1) - x(2:end) ; x(end)];
    DTx = @(x) ro*[x(1) ; x(2:end) - x(1 : end-1)];
    subg_F = @(x) ATA*x-ATb + (1/mu)*( DTx(Dx(x)) - DTx(Soft_Thresholding(Dx(x),mu)) );
    
    T = 10^8; %Time Horizon

    x_k = x_init;
    y_k = x_init;
    
    iter = 0;
    step = 1;
    L = L_f + alpha/mu;
    
    %First iteration f_k subg
    f_k = F(x_k);
    subg_k = subg_F(y_k);
    
    %
    fk_iter = zeros(T,1);
    fk_iter(iter+1) = f_k;
    for i=1:T-1
        if(f_k - f_opt<= epsilon)
            break;
        end
        
        new_x_k = y_k - (1/L)*subg_k;
        
        new_step = ( 1+sqrt(1+4*step^2) )/2 ;

        y_k = new_x_k  + ( (step-1)/new_step )*(new_x_k - x_k);
        
        x_k = new_x_k;
        step = new_step;
        
        %Next Iteration f_k subg
        f_k = F(x_k);
        subg_k = subg_F(y_k);
        iter = iter + 1;
        
        fk_iter(iter+1)=f_k;
    end
    fk_iter  = fk_iter(1:iter+1);
    x_est = x_k;
end