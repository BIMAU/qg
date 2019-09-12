function [x, l, k] = newtoncorrector(qg, par, ds, x, x0, l, l0, tol, varargin)

% Set some parameters
maxit = 20;
zeta = qg.get_par(20);

delta = 1e-6;

% Do the main iteration
for k = 1:maxit
    qg.set_par(par, l);
    fval = qg.rhs(x);

    qg.set_par(par, l + delta);
    dflval = (qg.rhs(x) - fval) / delta;
    qg.set_par(par, l);

    qg.jacob(x);

    z1 = -qg.solve(fval);
    z2 = qg.solve(dflval);

    rnp1 = zeta*(x-x0)'*(x-x0) + (1-zeta)*(l-l0)^2 - ds^2;
    dl = (-rnp1 - 2*zeta*(x-x0)'*z1) / (2*(1-zeta)*(l-l0) - 2*zeta*(x-x0)'*z2);
    dx = z1 - dl*z2;

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
