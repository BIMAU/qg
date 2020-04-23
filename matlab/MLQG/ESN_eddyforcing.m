rng(7)
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY X
% Load eddy forcing and coarse state data
% fname = 'eddyforcing_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend169_F0.5_Stir0_Rot1';
fname = 'eddyforcing_ff8_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend169_F0.5_Stir0_Rot1';
data  = load(['data/eddyforcing/', fname, '.mat']);
xaDim = size(data.xa,1);
efDim = size(data.eddyF,1);

% residual (eddy forcing) difference
data.rDiff = [zeros(efDim,1), data.eddyF(:,2:end) - data.eddyF(:,1:end-1)];

% Setup training data 
%  input  data U: averages and eddy forcings at time k*dt
%  output data Y: eddy forcings at time (k+1)*dt

N  = size(data.xa,2);  % total # samples
T  = 1600;             % training cutoff
Nt = T-1;              % # training samples
                       % Nt = 100;

assert(Nt <= T-1);

% create appropriate scalings for data
xaScaling = 3.0*max(abs(data.xa(:)));
efScaling = 2.0*max(abs(data.eddyF(:)));
rdScaling = 1.0*max(abs(data.rDiff(:)));

%scaleU = [xaScaling*ones(1,xaDim), efScaling*ones(1,efDim), rdScaling*ones(1,efDim)];
scaleU = [xaScaling*ones(1,xaDim); efScaling*ones(1,efDim)];
scaleU = scaleU(:)';
%scaleY = [efScaling*ones(1,efDim), rdScaling*ones(1,efDim)];
scaleY = [efScaling*ones(1,efDim)];

%U  = [data.xa(:,T-Nt:T-1); data.eddyF(:,T-Nt:T-1); data.rDiff(:,T-Nt:T-1)]';
U  = [data.xa(:,T-Nt:T-1); data.eddyF(:,T-Nt:T-1)]';
id = [(1:xaDim); (efDim+1:2*efDim)];
id = id(:)'; % alternate unknowns
U  = U(:,id);
Nu = size(U,2);

%Y  = [data.eddyF(:,T-Nt+1:T); data.rDiff(:,T-Nt+1:T)]';
Y  = [data.eddyF(:,T-Nt+1:T)]';
Ny = size(Y,2);

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;

tic; fprintf('create reservoir...\n');
createReservoir();
fprintf('create reservoir... done (%fs)\n', toc);

tic; fprintf('train reservoir...\n');
trainReservoir(U, Y);
fprintf('train reservoir... done (%fs)\n', toc);

%subplot(4,2,8)
%plot(X(2,:)); hold on
%plot(X(10,:)); hold on
%plot(X(end,:)); hold off
fprintf('begin prediction...\n');

% indices used for prediction (beyond T this is unseen data)
pMin = 10;
pPlus = 50;
pRange = T-pMin:T+pPlus;

% initial (known) eddy forcing
r  = data.eddyF(:,pRange(1)-1);

% initial (known) eddy forcing tangent
rd = data.rDiff(:,pRange(1)-1);

% reservoir state to use in autonomous mode
Rstate = X(end-pMin-1,:);

errnorm = [];
rdnorm  = [];
rdtnorm = [];
t = []; 
cl = lines(100);
cl = cl(iter,:);
iter = iter + 1;

%perturb data
%data.xa = data.xa.*(ones(xaDim, N) + 0.5*randn(xaDim, N));
for i = pRange
    xa = data.xa(:,i-1);
    %u  = [xa; r; rd];
    u  = [xa'; r'];
    u  = u(:);

    %Rstate(:) = update(Rstate, u' ./ scaleU, [r;rd]' ./ scaleY);
    Rstate(:) = update(Rstate, u' ./ scaleU, r' ./ scaleY);
    
    % predicted eddy forcing
    r0 = r;
    
    y  = tanh(W_out*[Rstate(:);u ./ scaleU']);
    y  = y .* scaleY'; % unscale
    r  = y(1:efDim);
    % rd = y(efDim+1:2*efDim);
    
    rd = r - r0;
    
    % r = r0 + rd;
    
    %-----------------------------------------------------------
    t  = [t, (data.times(i)-data.times(1)) / day];    
    fprintf('t = %d\n', round(t(end)));
    
    diff    = abs(r-data.eddyF(:,i));
    errnorm = [errnorm, norm(day*diff(:))];
    rdnorm  = [rdnorm, norm(day*rd(:))];
    rd_true = (data.eddyF(:,i) - data.eddyF(:,i-1));
    rdtnorm = [rdtnorm, norm(day*rd_true(:))];
    
    if ((i == pRange(end)) || (i == pRange(end)))
        subplot(4,2,1)
        plotQG(data.nxa, data.nya,1, day*data.xa(:,i), false)
        titleString = sprintf('Coarse vorticity day %d', round(t(end)));
        title(titleString);
        colorbar
        caxis([-0.2,0.2])

        subplot(4,2,2)
        plotQG(data.nxa, data.nya, 1, day*data.eddyF(:,i), false);
        colorbar
        caxis([-10,10])
        title('Eddy forcing')

        subplot(4,2,3)
        plotQG(data.nxa, data.nya, 1, day*r, false);
        colorbar
        caxis([-10,10])
        title('Predicted eddy forcing')

        subplot(4,2,4)
        plotQG(data.nxa, data.nya, 1, day*diff, false);
        colorbar    
        caxis([-10,10])
        if i < T
            col = 'k';
        else
            col = 'r';
        end
        
        title(sprintf('abs(diff), 2-norm = %3.2e', errnorm(end)),'color',col);
        
        subplot(4,2,5)
        plotQG(data.nxa, data.nya, 1, day*rd, false);
        colorbar
        %caxis([-10,10])
        title('rd (prediction)')


        subplot(4,2,6)
        plotQG(data.nxa, data.nya, 1, day*rd_true, false);
        colorbar
        %caxis([-10,10])
        title('rd (data)')
    end
    
    subplot(4,2,7)
    plot(t,errnorm,'.-','color',cl); hold on
    tsep = (data.times(T) - data.times(1)) / day;
    plot([tsep,tsep], ylim,'r-'); 

    subplot(4,2,8)
    plot(t,rdtnorm,'k.-'); hold on
    plot(t,rdnorm,'.-','color',cl); 
    tsep = (data.times(T) - data.times(1)) / day;
    plot([tsep,tsep], ylim,'r-');

    drawnow
end


function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY 
    % Reservoir parameters
    Nr       = 3000;
    noise    = 0.0;
    rhoMax   = 0.1;  % spectral radius
    entries  = 10;   % number of entries per row

    % create random matrix
    R = rand(Nr)-0.5;
    W = zeros(Nr);
    for i = 1:Nr
        [~,I] = sort(rand(Nr,1));
        I = I(1:10);
        W(i,I) = R(i,I);        
    end       

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
    Nu  = size(trainU,2);
    assert(dim == size(trainY, 1));
    % scale trainU, trainY
    %scaleU = 1.1*max(abs(trainU(:)));
    trainU = trainU ./ scaleU;
    %scaleY = 1.1*max(abs(trainY(:)));
    trainY = trainY ./ scaleY;

    % initialize activations X
    X = zeros(dim, Nr);
    
    subplot(4,2,2)
    plot(trainU(10,:));
    drawnow

    fprintf('#training fields: %d\n', dim);
    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end
    
    time = toc; fprintf('fitting W_out...\n')
    % compute W_out: W_out*X' = atanh(trainY)
    
    extX = [X, trainU];
    % Using a pseudo inverse
    % P = pinv(X);
    % W_out = (P*atanh(trainY))';
    
    % By solving the normal equations and including Tikhonov
    % regularization
    lambda  = 0.1;
    Xnormal = extX'*extX + lambda * speye(Nu+Nr);
    b       = extX'*atanh(trainY);
    W_out   = (Xnormal \ b)';
     
    fprintf('fitting W_out... done (%fs)\n', toc-time)
    % get training error
    predY = tanh(extX*W_out');

    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr
    pre = W*state' + W_in*u' + W_ofb*y';
    act = tanh(pre) + noise * (rand(Nr,1) - 0.5);
end