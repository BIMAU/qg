nx = 40;
ny = 40;
n = nx * ny * 2;
qg = QG(nx, ny);

% Set zeta
zeta = 1/n;
qg.set_par(20, zeta);

% List of parameters over which to perform a continuation
if ~exist('cont')
    cont = {};
    % Wind stress
    cont{1} = struct('par', 11, 'ds', 0.1, 'target', 1, 'maxit', 1000, 'from', [],...
                     'state', [], 'pars', [], 'history', []);
    % Asymmetry
    cont{2} = struct('par', 19, 'ds', 1, 'target', 1, 'maxit', 1, 'from', [],...
                     'state', [], 'pars', [], 'history', []);
    % Reynolds number
    cont{3} = struct('par', 5, 'ds', 1, 'target', 45, 'maxit', 1000, 'from', [],...
                     'state', [], 'pars', [], 'history', []);
    % Symmetry
    cont{4} = struct('par', 19, 'ds', -1, 'target', 0, 'maxit', 1000, 'from', [],...
                     'state', [], 'pars', [], 'history', []);
    % Reynolds number
    cont{5} = struct('par', 5, 'ds', -1, 'target', 45, 'maxit', 1000, 'from', [],...
                     'state', [], 'pars', [], 'history', []);
    % Reynolds number
    cont{6} = struct('par', 5, 'ds', 1, 'target', 45, 'maxit', 1000, 'from', 1,...
                     'state', [], 'pars', [], 'history', []);
end

x = zeros(n, 1);

x = newton(qg, x, 1e-10);

for i=1:length(cont)
    state = cont{i}.state;
    if ~isempty(state)
        for j=1:length(cont{i}.pars)
            qg.set_par(j, cont{i}.pars(j));
        end
        x = state;
        continue;
    else
        cont{i}.history = struct('psi_max', [], 'psi_min', [], 'par', []);
    end

    if ~isempty(cont{i}.from)
        x = cont{cont{i}.from}.state;
        for j=1:length(cont{cont{i}.from}.pars)
            qg.set_par(j, cont{cont{i}.from}.pars(j));
        end
    end

    par = cont{i}.par;
    ds = cont{i}.ds;
    target = cont{i}.target;
    maxit = cont{i}.maxit;

    % Get the initial tangent
    delta = 1e-6;
    l = qg.get_par(par);
    fval = qg.rhs(x);

    qg.set_par(par, l + delta);
    dl = (qg.rhs(x) - fval) / delta;
    qg.set_par(par, l);

    qg.jacob(x);
    dx = -qg.solve(dl);

    dl = 1;
    nrm = sqrt(zeta * dx'*dx + dl'*dl);
    dl = dl / nrm;
    dx = dx / nrm;

    dl0 = dl;
    dx0 = dx;

    % Perform the continuation
    for j=1:maxit
        l0 = l;
        x0 = x;

        l = l0 + ds * dl0;
        x = x0 + ds * dx0;

        [x2, l2] = newtoncorrector(qg, par, ds, x, x0, l, l0, 1e-10);

        dl = l2 - l0;
        l = l2;
        dx = x2 - x0;
        x = x2;

        psi_max = max(x(2:2:n));
        psi_min = min(x(2:2:n));
        fprintf('par = %2d, l = %f, max psi = %f, min psi = %f\n', ...
                par, l, psi_max, psi_min);
        cont{i}.history.psi_max = [cont{i}.history.psi_max, psi_max];
        cont{i}.history.psi_min = [cont{i}.history.psi_min, psi_min];
        cont{i}.history.par = [cont{i}.history.par, l];

        if (abs(dl) < 1e-10 || l >= target)
            % Converge onto the end point
            l = target;
            qg.set_par(par, l);
            x = newton(qg, x, 1e-10);

            psi_max = max(x(2:2:n));
            psi_min = min(x(2:2:n));
            fprintf('par = %2d, l = %f, max psi = %f, min psi = %f\n', ...
                    par, l, psi_max, psi_min);
            cont{i}.history.psi_max = [cont{i}.history.psi_max, psi_max];
            cont{i}.history.psi_min = [cont{i}.history.psi_min, psi_min];
            cont{i}.history.par = [cont{i}.history.par, l];
            break
        end

        dx0 = dx / ds;
        dl0 = dl / ds;
    end
    cont{i}.state = x;
    cont{i}.pars = [];
    for j=1:20
        cont{i}.pars = [cont{i}.pars, qg.get_par(j)];
    end
end

plot(cont{5}.history.par, cont{5}.history.psi_max,'.-');
hold on
plot(cont{6}.history.par, cont{6}.history.psi_max,'.-');
hold off