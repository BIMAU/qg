function root = newtonfw(f,df, x0, tol)

% Set some parameters
maxit = 100;
n = length(x0);

% Initialize output variable
x = x0;

% Do the main iteration
for k = 1:maxit
    fval = f(x);
    jac = df(x);
    dx = - jac\fval;
    x = x + dx;

    if norm(dx) < tol
        break
    end
end

root = x;

end
