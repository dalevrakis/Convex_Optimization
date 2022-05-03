function res = sgn(X)
    res = zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            if X(i,j) >= 0
                res(i,j) = 1;
            else
                res(i,j) = -1;
            end
        end
    end
end