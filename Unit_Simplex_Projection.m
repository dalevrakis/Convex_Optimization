function [x_proj] = Unit_Simplex_Projection(x)
    n=size(x,1);

    e = ones(n,1);

    eq = @(m) e'*max(x - m,0) - 1;
    
    % Bisection

    ub = max(x)+1; %upper bound
    lb = (sum(x)-1)/n - 1; %lower bound
    step = 0.1;

    MaxIter = 10^6;
    err = 10^-2;
    mu = 0;
    r = 1;

    % % Find good lower bound
    % while eq(lb) < 0
    %     lb = lb - step;
    % end
    % 
    % % Find good upper bound
    % while eq(ub) > 0
    %     ub = ub + step;
    % end

    iter = 0;
    for i = 1:MaxIter
        if(abs(r)<=err)
            break;
        end

        %Calc mid
        mu = (ub+lb)/2;
        r = eq(mu);
        if r > 0
            lb = mu;
        else
            ub = mu;
        end
        iter = iter + 1;
    end

    x_proj = max(x - mu,0);
end