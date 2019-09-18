clear all;
global nx ny nun n
nx= 64;
ny= 64;
nun=2;
n = nx*ny*nun;

idx1 = 1:nun:n;
idx2 = 2:nun:n;

global qg
qg = QG(nx, ny);

% enable wind forcing
qg.set_par(11,  1.0);
qg.set_par(5,  45.0);
qg.set_par(19,  0.01);


% mass matrix (remains cnst)
global B
B = qg.mass(n);

global history
history = struct('psi_max', [], 'psi_min', []);

% permutation matrix
bs = 4; % blocksize
P1 = blockpermutation(nx,ny,nun,bs);

% wavelet transform
H = haarmat(bs^2);
M = speye(n/(bs^2));
H = kron(M,H);

global Q nDetails nAverage

% permutation matrix
[P2, nDetails, nAverage] = separator(n, bs, 1);

% combined orthogonal transform
Q = P2*H*P1;

% time stepping
global t
t = 0;

global ZA ZD
ZA = [];
ZD = [];


% time stepping
dt = 0.01;
th = 1.0;

% initial solution
xold = zeros(n,1);

% number of time steps
Tend = 0.1;

useReservoir = false;
runTimeStepper(dt, th, Tend, xold)

global W W_in W_ofb W_out noise Nr scaleU scaleY
trainReservoir(ZA, ZD);





function [ ] = trainReservoir(trainU, trainY)
    global W W_in W_ofb W_out noise Nr scaleU scaleY

    Nr       = 100;
    noise    = 0.0;
    sparsity = .95;
    rhoMax   = 0.98;  % spectral radius

    W = rand(Nr)-0.5;
    W(rand(Nr) < sparsity) = 0;

    % Set spectral radius
    rho = eig(W);
    W   = W * rhoMax / max(abs(rho));

    % Create input weight matrix
    Nu = size(trainU,2);
    W_in = rand(Nr, Nu) * 2 - 1;

    % Create output feedback weight matrix
    Ny = size(trainY,2);
    W_ofb = rand(Nr, Ny) * 2 - 1;

    dim = size(trainU,1);

    % scale trainU, trainY
    scaleU = 0.9*max(abs(trainU(:)));
    trainU = trainU / scaleU;
    scaleY = 0.9*max(abs(trainY(:)));
    trainY = trainY / scaleY;

    % initialize activations X
    X = zeros(dim, Nr);

    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end

    % Extend X with raw input columns
    % extX = [X, trainU];
    extX = X;

    state = X(end,:);

    % Create pseudo inverse
    P = pinv(extX);

    % compute W_out
    W_out = (P*atanh(trainY))';
    predY = tanh(extX*W_out');
    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

end

% overload rhs function with ML component
function [out] = F(x)
    global W_out scaleU Q nDetails nAverage state useReservoir qg

    if useReservoir
        z     = Q*x;
        z     = scaleU * z(1:nAverage)';
        state = update(state, z, z);

        y = tanh(W_out*z');
        y = y / scaleU;
        z = [z';y'];
        x = Q'*z;
    end
    % TODO
    out = qg.rhs(x);
end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr
    pre = W*state' + W_in*u' + W_ofb*y';
    act = tanh(pre) + noise * (rand(Nr,1) - 0.5);
end

function [] = runTimeStepper(dt, th, Tend, x0)
    global B qg history ZA ZD t Q nAverage nDetails
    global nx ny nun n

    % initial state
    xold = x0;

    % initial rhs
    global useReservoir
    useReservoir = false;
    rhsold = F(xold);

    % scaling for time dependent component of Jacobian
    sig  = 1.0/(dt*th);

    % optimal number of Newton steps
    Nopt = 3.5;

    % max number of Newton steps
    nsteps = 10;

    while t < Tend
        t = t + dt;
        x = xold;
        for k = 1:nsteps
            rhs = B*(x-xold)/(dt*th) + F(x) + (1-th)/th * rhsold;
            qg.jacob(x, sig);
            dx = qg.solve(-rhs);
            x = x + dx;
            if norm(dx) < 1e-6
                fprintf('Newton corrector converged in %d steps\n', k);
                break;
            end
        end
        if k == nsteps
            fprintf('Warning: Newton failed, residual =  %f\n', ...
                    norm(dx));
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

        psi_max = max(x(2:2:end));
        psi_min = min(x(2:2:end));

        history.psi_max = [history.psi_max, psi_max];
        history.psi_min = [history.psi_min, psi_min];

        fprintf('time step size: %f\n', dt);
        fprintf('max(psi) = %e\n', psi_max);
        fprintf('time t = %f\n', t);

        % store data for next time step
        xold   = x;
        rhsold = F(x);

        z  = Q*x;
        za = z; za(nAverage+1:end) = 0;
        zd = z; zd(1:nAverage) = 0;
        ZA = [ZA; z(1:nAverage)'];
        ZD = [ZD; z(nAverage+1:end)'];

        xa = Q'*za;
        xd = Q'*zd;

        figure(1); plotQG(nx,ny,2,xd); drawnow;
        figure(2); plotQG(nx,ny,2,xa); drawnow;
        figure(3); plotQG(nx,ny,2,x);  drawnow;
    end

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