function root = newton(qg, x0, tol, varargin)

% Set some parameters
maxit = 1000;
n = length(x0);

% Initialize output variable
x = x0;

% Do the main iteration
for k = 1:maxit
    fval = qg.rhs(x);
    qg.jacob(x);
    dx = -qg.solve(fval);

    x = x + dx;

    if norm(dx) < tol
        break
    end
end

root = x;

end
