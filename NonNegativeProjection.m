function Px = NonNegativeProjection(x)
    
    Px = zeros(size(x));
    
    for i = 1:size(x,1)
        if x(i)<0
            Px(i) = 0;
        else
            Px(i) = x(i);
        end
    end
    
end