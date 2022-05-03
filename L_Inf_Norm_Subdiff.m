function subdiff = L_Inf_Norm_Subdiff(A,b,x)
    n = size(x,1);
    m = size(A,1);
%     subdiff = zeros(n,1);
    res = A*x-b;
    if x == zeros(n,1)
        subdiff = (1/m)*A'*ones(m,1);
    else
        I_x = find(abs(res)== max(abs(res)));
        s_I = size(I_x,1);
        subdiff = zeros(n,s_I);
        
        for i=1:s_I
            I_ind = I_x(i);
            subdiff(:,i) = sgn( A(I_ind,:)*x-b(I_ind) )*A(I_ind,:)';
        end
    end
end