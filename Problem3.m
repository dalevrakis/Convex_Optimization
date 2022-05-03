n = 4;

e = ones(n,1);

c = rand(n,1);

cvx_begin
    variables x(n)
    minimize (c'*x);
    subject to
        e'*x-1==0;
        x>=0;
cvx_end