function UB = Upper_Bound(iter,L_f,x_0,f_opt)
    UB = zeros(iter,1);
    for i=1:iter
        UB(i) = (L_f*norm(x_0-f_opt,2))/sqrt(i+1);
    end
end