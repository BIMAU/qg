
nx= 64;
ny= 64;
n = nx*ny*2;

idx1 = 1:2:n;
idx2 = 2:2:n;

qg = QG(nx, ny);

% enable wind forcing
qg.set_par(11,  1.0);
qg.set_par(5,  45.0);
qg.set_par(19,  0.01);

% time stepping
dt = 0.01;
th = 1.0;

% scaling for time dependent component of Jacobian
sig  = 1.0/(dt*th);

% mass matrix (remains cnst)
B = qg.mass(n);

% initial solution
xold = zeros(n,1);

% initial rhs
rhsold = qg.rhs(xold);

% number of time steps
Tend = 1000;

% max number of Newton steps
nsteps = 20;

% optimal number of Newton steps
Nopt = 3.5;

history = struct('psi_max', [], 'psi_min', []);

% time stepping
t = 0;
while t < Tend
    t = t + dt;
    x = xold;
    for k = 1:nsteps
        rhs = B*(x-xold)/(dt*th) + qg.rhs(x) + (1-th)/th * rhsold;
        qg.jacob(x, sig);
        dx = qg.solve(-rhs);
        x = x + dx;
        if norm(dx) < 1e-6
            fprintf('Newton corrector converged in %d steps\n', k);
            break;
        end
    end
    
    % step size control
    factor = Nopt / k;
    factor = max(0.5, factor);
    factor = min(2.0, factor);    
    dt   =  dt * factor;
         
    if dt > 1000
        fprintf('Equilibrium reached.\n', k);
        break;
    end
    sig  = 1.0/(dt*th);
    
    psi_max = max(x(2:2:n));
    psi_min = min(x(2:2:n));
    
    history.psi_max = [history.psi_max, psi_max];
    history.psi_min = [history.psi_min, psi_min];
    
    fprintf('time step size: %f\n', dt);
    fprintf('max(psi) = %e\n', psi_max);
    fprintf('time t = %f\n', t);
    
    % store data for next time step
    xold   = x;
    rhsold = qg.rhs(x);       
    
    plotQG(nx,ny,2,x); drawnow
end















% $$$ % test analytical jacobian with numerical jacobian
% $$$ J = -qg.jacobian(xold, 0.0);
% $$$ 
% $$$ pert = 1e-6;
% $$$ err = zeros(n,1);
% $$$ for i = 1:n
% $$$     xold(i) = xold(i) + pert;    
% $$$     err(i)  = max(abs(( qg.rhs(xold) - rhsold ) / pert  - J(:,i) ));    
% $$$     xold(i) = xold(i) - pert;   
% $$$ end
% $$$ 
% $$$ plot(err);hold on
% $$$ norm(err)
% $$$ 
% $$$ % test jacobian with time dependent contribution with numerical jacobian

% $$$ G = @(x) -B*(x-xold)/(dt*th) + qg.rhs(x) + (1-th)/th * rhsold
% $$$ err    = zeros(n,1);
% $$$ rhsold = G(xold);
% $$$ J = -qg.jacobian(xold, sig);
% $$$ for i = 1:n
% $$$     xold(i) = xold(i) + pert;    
% $$$     err(i)  = max(abs(( G(xold) - rhsold ) / pert  - J(:,i) ));    
% $$$     xold(i) = xold(i) - pert;   
% $$$ end
% $$$ plot(err); hold off
% $$$ norm(err)