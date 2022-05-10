function T_l = Soft_Thresholding(x, lambda)
    T_l = max(abs(x)-lambda,0).*sgn(x);
end