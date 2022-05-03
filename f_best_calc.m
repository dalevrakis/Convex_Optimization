function f_best = f_best_calc(f_k)
    s = size(f_k,1);
    f_best = zeros(s,1);
   
    f_min = f_k(1);
    f_best(1) = f_min;
    for i=2:s
        if f_k(i) < f_min
            f_min = f_k(i);
        end
        f_best(i) = f_min;
    end
end