function [x, l] = newtoncorrector(qg, par, ds, x, x0, l, l0, tol, varargin)

% Set some parameters
maxit = 20;
zeta = qg.get_par(20);

delta = 1e-6;

% Do the main iteration
for k = 1:maxit
    % Set the parameter value and compute F (RHS of 2.2.9)
    qg.set_par(par, l);
    fval = qg.rhs(x);

    % Compute F_mu (bottom part of the RHS of 2.2.9)
    qg.set_par(par, l + delta);
    dflval = (qg.rhs(x) - fval) / delta;
    qg.set_par(par, l);

    % Compute the jacobian at x
    qg.jacob(x);

    % Solve twice with F_x (2.2.9)
    z1 = -qg.solve(fval);
    z2 = qg.solve(dflval);

    % Compute r (2.2.8)
    rnp1 = zeta*(x-x0)'*(x-x0) + (1-zeta)*(l-l0)^2 - ds^2;

    % Compute dl (2.2.13)
    dl = (-rnp1 - 2*zeta*(x-x0)'*z1) / (2*(1-zeta)*(l-l0) - 2*zeta*(x-x0)'*z2);

    % Compute dx (2.2.12)
    dx = z1 - dl*z2;

    % Compute a new x and l (2.2.10 - 2.2.11)
    x = x + dx;
    l = l + dl;

    % fprintf('%f, %e, %e\n', l, dl, norm(dx));
    if norm(dx) < tol
        fprintf('Newton corrector converged in %d steps\n', k);
        return
    end
end
fprintf('No convergence achieved by Newton corrector\n');

end
