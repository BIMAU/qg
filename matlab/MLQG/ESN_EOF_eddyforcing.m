rng(7)
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY X alpha

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;

% Load eddy forcing and coarse state data
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend151_F0.5_Stir0_Rot1';
% data  = load(['data/eddyforcing/', fname, '.mat']);
M     = size(data.xa,1);
assert(M == size(data.eddyF,1));

nxa = data.nxa;
nya = data.nya;
nun = 2;
assert(M == nxa*nya*nun);

dimx = 64;
dimr = 2;
% [Vx, ~] = pca(data.xa, dimx);
% [Vr, ~] = pca(data.eddyF, dimr);


N  = size(data.xa,2);  % total # samples/experiments
T  = 3000;             % training cutoff
Nt = T-1;              % # training samples
                       % Nt = 100;
                       % Nt = 500;

assert(Nt <= T-1);

xbar = Vx' * data.xa(:,T-Nt:T-1);
rbar = Vr' * data.eddyF(:, T-Nt:T-1);
rbrp = Vr' * data.eddyF(:, T-Nt+1:T);
mxx  = 3.0*max(abs(xbar(:)));  % why does this work so well??
mxr  = 3.0*max(abs(rbar(:)));

% U = [ones(1,Nt); xbar; rbar]';
U = [xbar; rbar]';
% scaleU = [1, mxx*ones(1, dimx), mxr*ones(1, dimr)];
scaleU = [mxx*ones(1, dimx), mxr*ones(1, dimr)];
% Nu = 1 + dimx + dimr;
Nu = dimx + dimr;

Y = rbrp';
scaleY = mxr*ones(1, dimr);
Ny = dimr;

alpha = 0.5;
tic; fprintf('create reservoir...\n');
createReservoir();
fprintf('create reservoir... done (%fs)\n', toc);

tic; fprintf('train reservoir...\n');
trainReservoir(U, Y);
fprintf('train reservoir... done (%fs)\n', toc);

fprintf('begin prediction...\n');

% indices used for prediction (beyond T this is unseen data)
pMin   = 10;
pPlus  = 30;
pRange = T-pMin:T+pPlus;

% initial (known) eddy forcing
r  = Vr' * data.eddyF(:,pRange(1)-1);

% reservoir state to use in autonomous mode
Rstate = X(end-pMin-1,:);

errnorm = [];
rdnorm  = [];
rdtnorm = [];
t  = []; 
cl = lines(100);
cl = cl(iter,:);
iter = iter + 1;

for i = pRange
    xa = Vx' * data.xa(:,i-1);

    %u  = [1; xa; r];
    u  = [xa; r];

    Rstate(:) = update(Rstate, u' ./ scaleU, r' ./ scaleY);
    
    % predicted eddy forcing
    
    y  = tanh(W_out*[Rstate(:); u ./ scaleU']);
    % y  = tanh(W_out*Rstate(:));
    y  = y .* scaleY'; % unscale
    r  = y(1:dimr);    
    
    
    % postprocessing ---------------------------------
    diff    = abs(r-(Vr'*data.eddyF(:,i)));
    errnorm = [errnorm, norm(diff(:))];
    t       = [t, (data.times(i)-data.times(1)) / day];    
    
    subplot(2,2,1)
    plotQG(nxa, nya, 1, day*Vr*(Vr'*data.eddyF(:,i)), false)
    title('projected eddy forcing')

    subplot(2,2,2)
    %plotQG(nxa, nya, 1, day*V*(V'*data.eddyF(:,i)), false)
    %caxis([-10,10]);
    plot((Vr'*data.eddyF(:,i)),'.-'); hold on;
    plot(r,'.-'); hold off
    legend('true', 'ML')
    title('eddy forcing coefficients')
    
    subplot(2,2,3)
    plotQG(nxa, nya, 1, day*Vr*r, false)    
    title('predicted eddy forcing')
    
    subplot(2,2,4)
    plot(t,errnorm,'.-','color',cl); hold on
    tsep = (data.times(T) - data.times(1)) / day;
    plot([tsep,tsep], ylim,'r-'); 
    ylim([0,5000]);
    
    drawnow
end


function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY 
    % Reservoir parameters
    Nr       = 1000;
    noise    = 0.0;
    rhoMax   = 0.90;  % spectral radius
                     %    sparsity = 0.;
    entries  = 10;   % number of entries per row

    % create random matrix
    R = rand(Nr)-0.5;
    W = zeros(Nr);
    for i = 1:Nr
        [~,I] = sort(rand(Nr,1));
        I = I(1:10);
        W(i,I) = R(i,I);        
    end       
    
    % W = rand(Nr)-0.5;
    % W(rand(Nr) < sparsity) = 0;
    rho = eig(W);
    W   = W * rhoMax / max(abs(rho));
    
    % Set spectral radius
    % rho = eigs(W, 1);
    % W   = W * rhoMax / abs(rho);

    % Create input weight matrix
    % D = [(1:Nr)', ceil(Nu*rand(Nr,1)), (rand(Nr,1)*2-1)];
    % W_in = spconvert(D);
    W_in = (rand(Nr, Nu) * 2 - 1);
    
    W_in(1:Nr/2, Nu/2+1:end) = 0;
    W_in(Nr/2+1:end, 1:Nu/2) = 0;

    % Create output feedback weight matrix
    % D = [(1:Nr)', ceil(Ny*rand(Nr,1)), (rand(Nr,1)*2-1)];
    % W_ofb = spconvert(D);
    W_ofb = (rand(Nr, Ny) * 2 - 1);    

end


function [ ] = trainReservoir(trainU, trainY)
    global W W_in W_ofb W_out noise Nr scaleU scaleY X
    
    dim = size(trainU,1);
    Nu  = size(trainU,2);
    assert(dim == size(trainY, 1));

    % scale trainU, trainY
    trainU = trainU ./ scaleU;
    trainY = trainY ./ scaleY;

    % initialize activations X
    X = zeros(dim, Nr);
    
    fprintf('#training fields: %d\n', dim);
    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end
    
    time = toc; fprintf('fitting W_out...\n')
    % compute W_out: W_out*X' = atanh(trainY)
    
    % extX = X;
    extX = [X, trainU];

    % Using a pseudo inverse
    % P = pinv(X);
    % W_out = (P*atanh(trainY))';
    
    % By solving the normal equations and including Tikhonov
    % regularization
    lambda  = 0.01;
    Xnormal = extX'*extX + lambda * speye(Nu+Nr);
    % Xnormal = extX'*extX + lambda * speye(Nr);
    
    b       = extX'*atanh(trainY);
    W_out   = (Xnormal \ b)';
     
    fprintf('fitting W_out... done (%fs)\n', toc-time)

    % get training error
    predY = tanh(extX*W_out');

    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr alpha
    pre = W*state' + W_in*u' + W_ofb*y';
    act = alpha * tanh(pre) + (1-alpha) * state' + noise * (rand(Nr,1) - 0.5);
end