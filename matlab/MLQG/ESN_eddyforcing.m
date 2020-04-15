
rng(7);
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY X
% Load eddy forcing and coarse state data
fname = 'eddyforcing_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
data  = load(['data/eddyforcing/', fname, '.mat']);

% Setup training data 
%  input  data U: averages and eddy forcings at time k*dt
%  output data Y: eddy forcings at time (k+1)*dt

N  = size(data.xa,2);  % total # samples
T  = 220;              % training cutoff
Nt = T-1;              % # training samples
Nt = 100;
assert(Nt <= T-1);

% create appropriate scalings for data
xaScaling = 2.0*max(abs(data.xa(:)));
efScaling = 2.0*max(abs(data.eddyF(:)));
xaDim = size(data.xa,1);
efDim = size(data.eddyF,1);

scaleU = [xaScaling*ones(1,xaDim), efScaling*ones(1,efDim)];
scaleY = [efScaling*ones(1,efDim)];

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
tic
trainReservoir(U, Y);
fprintf('train reservoir... done (%fs)\n', toc);
subplot(4,2,8)
plot(X(2,:)); hold on
plot(X(10,:)); hold on
plot(X(end,:)); hold off
fprintf('begin prediction...\n');
% indices used for prediction (beyond T this is unseen data)
pMin = 10;
pPlus = 5;
pRange = T-pMin:T+pPlus;

% initial (known) eddy forcing
r  = data.eddyF(:,pRange(1)-1);

% reservoir state to use in autonomous mode
Rstate = X(end-pMin,:);

errnorm = [];
t = []; 
cl = lines(100);
cl = cl(iter,:);
iter = iter + 1;
for i = pRange
    xa = data.xa(:,i-1);
    u  = [xa; r];

    Rstate(:) = update(Rstate, u' ./ scaleU, r' ./ scaleY);
    
    % predicted eddy forcing
    r0 = r;
    r = tanh(W_out*Rstate(:));
    r = r .* scaleY'; % unscale
    
    %-----------------------------------------------------------
    t  = [t, (data.times(i)-data.times(1)) / day];    
    fprintf('t = %f\n', t(end));
    
    figure(1)
    subplot(4,2,1)
    plotQG(data.nxa, data.nya,1, day*data.xa(:,i), false)
    titleString = sprintf('Coarse vorticity day %f', t(end));
    title(titleString);
    colorbar
    caxis([-0.2,0.2])

    subplot(4,2,2)
    plotQG(data.nxa, data.nya, 1, day*data.eddyF(:,i), false);
    colorbar
    caxis([-10,10])
    title('Eddy forcing for coarse model')

    subplot(4,2,3)
    plotQG(data.nxa, data.nya, 1, day*r, false);
    colorbar
    caxis([-10,10])
    title('Predicted eddy forcing for coarse model')

    diff = abs(r-data.eddyF(:,i));
    subplot(4,2,4)
    plotQG(data.nxa, data.nya, 1, day*diff, false);
    colorbar    
    caxis([-10,10])
    if i < T
        col = 'k';
    else
        col = 'r';
    end
    
    errnorm = [errnorm, norm(day*diff(:))];
    title(sprintf('abs(diff), 2-norm = %3.2e', errnorm(end)),'color',col);
    
    drdt = (r - r0);% / (data.times(i)-data.times(i-1));
    subplot(4,2,5)
    plotQG(data.nxa, data.nya, 1, day*drdt, false);
    colorbar
    %caxis([-10,10])
    title('drdt (prediction)')

    drdt_true = (data.eddyF(:,i) - data.eddyF(:,i-1));% / (data.times(i)-data.times(i-1));
    subplot(4,2,6)
    plotQG(data.nxa, data.nya, 1, day*drdt_true, false);
    colorbar
    %caxis([-10,10])
    title('drdt (data)')
    
    subplot(4,2,7)
    plot(t,errnorm,'.-','color',cl); hold on
    tsep = (data.times(T) - data.times(1)) / day;
    plot([tsep,tsep], ylim,'r-');
    
    drawnow
end


function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY 
    % Reservoir parameters
    Nr       = 300;
    noise    = 0.;
    sparsity = 0.92;
    rhoMax   = 1.0;  % spectral radius

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
    global W W_in W_ofb W_out noise Nr scaleU scaleY X
    
    dim = size(trainU,1);

    % scale trainU, trainY
    %scaleU = 1.1*max(abs(trainU(:)));
    trainU = trainU ./ scaleU;
    %scaleY = 1.1*max(abs(trainY(:)));
    trainY = trainY ./ scaleY;

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
    lambda  = 1; 
    Xnormal = X'*X + lambda * speye(Nr);
    b       = X'*atanh(trainY);
    W_out   = (Xnormal \ b)';
     
    % get training error
    predY = tanh(X*W_out');
    fprintf('#training fields: %d\n', dim);
    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr
    pre = W*state' + W_in*u' + W_ofb*y';
    act = tanh(pre) + noise * (rand(Nr,1) - 0.5);
end