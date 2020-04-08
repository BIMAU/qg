clear all;
rng(7);
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY Rstate
% Load eddy forcing and coarse state data
fname = 'eddyforcing_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
data  = load(['data/eddyforcing/', fname, '.mat']);

% Setup training data 
%  input  data U: averages and eddy forcings at time k*dt
%  output data Y: eddy forcings at time (k+1)*dt

N  = size(data.xa,2);  % total # samples
T  = 100;              % training cutoff
Nt = 50;               % # training samples

U  = [data.xa(:,T-Nt:T-1); data.eddyF(:,T-Nt:T-1)]';
Nu = size(U,2);
Y  = [data.eddyF(:,T-Nt+1:T)'];
Ny = size(Y,2);

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;

fprintf('create reservoir...\n');
createReservoir();
fprintf('create reservoir... done\n');

fprintf('train reservoir...\n');
trainReservoir(U, Y);
fprintf('train reservoir... done\n');

fprintf('begin prediction...\n');
% indices used for prediction (beyond T this is unseen data)
pRange = T-10:T+10;

% initial (known) eddy forcing
r  = data.eddyF(:,pRange(1)-1);
for i = pRange
    xa = data.xa(:,i-1);
    u  = [xa; r];

    Rstate(:) = update(Rstate, u' / scaleU , r' / scaleY);
    
    % predicted eddy forcing
    r0 = r;
    r = tanh(W_out*Rstate(:));
    r = r * scaleY; % unscale
    
    %-----------------------------------------------------------
    t  = data.times(i);    
    subplot(3,2,1)
    plotQG(data.nxa, data.nya,1, day*data.xa(:,i), false)
    titleString = sprintf('Coarse vorticity day %f', (t-data.times(1)) / day);
    title(titleString);
    colorbar
    caxis([-0.2,0.2])

    subplot(3,2,2)
    plotQG(data.nxa, data.nya, 1, day*data.eddyF(:,i), false);
    colorbar
    caxis([-10,10])
    title('Eddy forcing for coarse model')

    subplot(3,2,3)
    plotQG(data.nxa, data.nya, 1, day*r, false);
    colorbar
    caxis([-10,10])
    title('Predicted eddy forcing for coarse model')

    diff = abs(r-data.eddyF(:,i));
    subplot(3,2,4)
    plotQG(data.nxa, data.nya, 1, day*diff, false);
    colorbar    
    caxis([-10,10])
    if i < T
        col = 'k';
    else
        col = 'r';
    end
    title(sprintf('abs(diff), 2-norm = %3.2f', norm(day*diff(:))),'color',col);
    
    drdt = (r - r0) / (data.times(i)-data.times(i-1));
    subplot(3,2,5)
    plotQG(data.nxa, data.nya, 1, day*drdt, false);
    colorbar
    %caxis([-10,10])
    title('drdt (prediction)')

    drdt_true = (data.eddyF(:,i) - data.eddyF(:,i-1)) / (data.times(i)-data.times(i-1));
    subplot(3,2,6)
    plotQG(data.nxa, data.nya, 1, day*drdt_true, false);
    colorbar
    %caxis([-10,10])
    title('drdt (data)')

    
    drawnow
    %exportfig('out.eps',10,[18,18]);

end


function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY Rstate
    % Reservoir parameters
    Nr       = 1000;
    noise    = 0.0;
    sparsity = 0.9;
    rhoMax   = 0.9;  % spectral radius

    % create random matrix
    W = rand(Nr)-0.5;

    % set sparsity
    W(rand(Nr) < sparsity) = 0;

    % Set spectral radius
    rho = eig(W);
    W   = W * rhoMax / max(abs(rho));
    W   = sparse(W);

    % Create input weight matrix
    W_in = rand(Nr, Nu) * 2 - 1;

    % Create output feedback weight matrix
    W_ofb = rand(Nr, Ny) * 2 - 1;
end

function [ ] = trainReservoir(trainU, trainY)
    global W W_in W_ofb W_out noise Nr scaleU scaleY Rstate
    
    dim = size(trainU,1);

    % scale trainU, trainY
    scaleU = 1.1*max(abs(trainU(:)));
    trainU = trainU / scaleU;
    scaleY = 1.1*max(abs(trainY(:)));
    trainY = trainY / scaleY;

    % initialize activations X
    X = zeros(dim, Nr);

    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end

    % compute W_out: W_out*X' = atanh(trainY)
    
    % Using a pseudo inverse
    % P = pinv(X);
    % W_out = (P*atanh(trainY))';
    
    % By solving the normal equations and including Tikhonov
    % regularization
    lambda  = 100; 
    Xnormal = X'*X + lambda * speye(Nr);
    b       = X'*atanh(trainY);
    W_out   = (Xnormal \ b)';
     
    % get training error
    predY = tanh(X*W_out');
    fprintf('#training fields: %d\n', dim);
    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

    % save last reservoir state to use in autonomous mode
    Rstate = X(end,:);
end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr
    pre = W*state' + W_in*u' + W_ofb*y';
    act = tanh(pre) + noise * (rand(Nr,1) - 0.5);
end