rng(77)
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY X alpha x1max x2max nxa nya 

global extendX outputActivation reservoirStateInit

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;

% Load eddy forcing and coarse state data
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend151_F0.5_Stir0_Rot1';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart151_Tend179_F0.5_Stir0_Rot1';
% data  = load(['data/eddyforcing/', fname, '.mat']);
M = size(data.xa, 1);
S = size(data.xa, 2);

% block ordering and scaling
x1    = data.xa(1:2:end,:);
x1max = max(abs(x1(:)));
x2    = data.xa(2:2:end,:);
x2max = max(abs(x2(:)));
data.blockxa = [x1 / x1max; x2 / x2max];

assert(M == size(data.eddyF,1));

nxa = data.nxa;
nya = data.nya;
nun = 2;
assert(M == nxa*nya*nun);

sep  = 0;
dimx = 128-sep;
dimr = 2;

extendX = true;
outputActivation = false;
reservoirStateInit = 'rand';  % 'zero', 'rand'

% [V, ~,var] = pca(data.blockxa);
% [Vr, ~] = pca(data.eddyF, dimr);
Vx  = V(:,sep+1:sep+dimx);
Va  = V(:,1:sep);
N   = size(data.xa,2);  % total # samples/experiments
T   = round(0.65*S);     % training cutoff
Nt  = T-1;              % # training samples
                        % Nt = 100;

% Nt = 2000;

assert(Nt <= T-1);

proj  = [Va,Vx]' * data.blockxa;
mxx   = 4. * max(abs(proj(:)));
abar  = Va' * data.blockxa(:,T-Nt:T-1);
abarp = Va' * data.blockxa(:,T-Nt+1:T);
xbar  = Vx' * data.blockxa(:,T-Nt:T-1);
xbarp = Vx' * data.blockxa(:,T-Nt+1:T);
% rbar  = Vr' * data.eddyF(:, T-Nt:T-1);
% rbarp = Vr' * data.eddyF(:, T-Nt+1:T);
mxa = mxx;
mxx = mxx;
% mxr  = 1.2*max(abs(rbar(:)));

% U = [ones(1,Nt); xbar; rbar]';
% U = [xbar; rbar]';
U = [abar; xbar]';
% U = [abar]';
% scaleU = [1, mxx*ones(1, dimx), mxr*ones(1, dimr)];
% scaleU = [mxx*ones(1, dimx), mxr*ones(1, dimr)];
scaleU = [mxa*ones(1,sep), mxx*ones(1, dimx)];
% scaleU = [mxa*ones(1,sep)];

% Nu = 1 + dimx + dimr;
% Nu = dimx + dimr;
Nu = sep+dimx;
% Nu = sep;

% Y = rbarp';
Y = [abarp; xbarp]';
% Y = xbarp';
% scaleY = mxr*ones(1, dimr);
scaleY = [mxa*ones(1,sep), mxx*ones(1, dimx)];
% scaleY = [mxx*ones(1, dimx)];
% Ny = dimr;
Ny = sep+dimx;
% Ny = dimx;

alpha = 0.7;
beta  = 1.0;

tic; fprintf('create reservoir...\n');
createReservoir();
fprintf('create reservoir... done (%fs)\n', toc);

tic; fprintf('train reservoir...\n');
trainReservoir(U, Y);
figure(2)
colormap(gray)
subplot(2,1,1)
imagesc(X(:,Nr/2:Nr/2+1+10)');
colorbar
caxis([-1,1])
title('Reservoir')
subplot(2,1,2)
imagesc(X(:,:)');
colorbar
caxis([-1,1])
fprintf('train reservoir... done (%fs)\n', toc);

fprintf('begin prediction...\n');
% indices used for prediction (beyond T this is unseen data)
pMin   = 0;
pPlus  = S-T;
pRange = T-pMin:T+pPlus;

% initial (known) eddy forcing
% r  = Vr' * data.eddyF(:,pRange(1)-1);
xa = Va' * data.blockxa(:,pRange(1)-1);
x  = Vx' * data.blockxa(:,pRange(1)-1);
xp = Vx' * data.blockxa(:,pRange(1));

% reservoir state to use in autonomous mode
Rstate = X(end-pMin-1,:);
Xpred  = zeros(numel(pRange), Nr);

errnorm = [];
conorm1 = [];
conorm2 = [];
rdnorm  = [];
rdtnorm = [];
t  = [];
cl = lines(100);
cl = cl(iter,:);
iter = iter + 1;
dt = data.times(2)-data.times(1);

figure(1)
for k = 1:numel(pRange)
    i   = pRange(k);
    xat = beta * Va' * data.blockxa(:,i-1) + (1-beta)*xa;
    % u  = [1; xa; r];
    % u  = [xa; r];
    u = [xat; x];
    % u = xa;
    y = [xat; x];
    % y = x;

    Rstate(:) = update(Rstate, u' ./ scaleU, y' ./ scaleY);
    Xpred(k,:) = Rstate;
    
    % prediction    
    if extendX
        y  = tanh(W_out*[Rstate(:); u ./ scaleU']);
    else
        y  = tanh(W_out*Rstate(:));
    end
    
    y = y .* scaleY';  % unscale

    x  = y(sep+1:sep+dimx);
    xa = y(1:sep);
        
    % postprocessing ---------------------------------
    co1     = ([Va,Vx]'*data.blockxa(:,i));
    co2     = [xat; x];
    diff    = abs(co2-co1);
    errnorm = [errnorm, norm(diff(:))];
    t       = [t, (data.times(i)-data.times(1)) / day];    
    conorm1 = [conorm1, norm(co1(:))];
    conorm2 = [conorm2, norm(co2(:))];
    
    if ((mod(i,1000) == 0) || (i == pRange(end)))
        subplot(3,2,1)
        plotBlockQG(day*[Va,Vx]*([Va,Vx]'*data.blockxa(:,i)));
        title('Model, vorticity (1/day)')
        caxis([-0.5,0.5])

        subplot(3,2,2)
        plot(([Va,Vx]'*data.blockxa(:,i)),'.-'); hold on;
        plot([xat;x],'.-'); hold off
        legend('Model', 'Reservoir')
        title('EOF coefficients')
        xlim([1,min(100,Ny)])
        
        subplot(3,2,3)
        plotBlockQG(day*[Va,Vx]*[xat;x])
        title('Reservoir prediction')
        caxis([-0.5,0.5])
        
        subplot(3,2,4)    
        plot(t,conorm2.^2,'.-','color',cl); hold on
        plot(t,conorm1.^2,'k-');
        title('Energy')
        
        subplot(3,2,5)
        imagesc(Xpred');
        colorbar
        caxis([-1,1])
        title('Reservoir');
        
        subplot(3,2,6)
        hd(iter) = plot(t,errnorm,'.-','color',cl); hold on
        tsep = (data.times(T) - data.times(1)) / day;
        % semilogy([tsep,tsep], [1e-1,1e2],'r-');
        title('Error norm')
        drawnow         
    end
end

function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY 
    % Reservoir parameters
    Nr     = 4000;
    noise  = 0.0;
    rhoMax = 0.95;  % spectral radius
                    % sparsity = 0.
    
    entriesPerRow = 10;
    
    D = [];
    for i = 1:entriesPerRow
        D = [D; [(1:Nr)', ceil(Nr*rand(Nr,1)), (rand(Nr,1)-0.5)] ];
    end
    W = spconvert(D);
    
    % W = rand(Nr)-0.5;
    % W(rand(Nr) < 0.98) = 0;

    % Set spectral radius
    opts.maxit=1000;
    rho = eigs(W, 2, 'lm', opts);
    mrho = max(abs(rho));
    if isnan(mrho)
        throw();
    end
    fprintf('using max eig: %f\n', mrho);
    W   = W * rhoMax / mrho;

    % Create input weight matrix
    D = [(1:Nr)', ceil(Nu*rand(Nr,1)), (rand(Nr,1)*2-1)];
    W_in = sparse(D(:,1),D(:,2),D(:,3),Nr,Nu);
    
    % W_in(rand(Nr,1)<0.5,:) = 0.0;
    % W_in = speye(Nr,Nu);
    % W_in = (rand(Nr, Nu) * 2 - 1);
    
    % W_in(1:Nr/2, Nu/2+1:end) = 0;
    % W_in(Nr/2+1:end, 1:Nu/2) = 0;

    % Create output feedback weight matrix
    % D = [(1:Nr)', ceil(Ny*rand(Nr,1)), (rand(Nr,1)*2-1)];
    % W_ofb = spconvert(D);
    W_ofb = 0.1*(rand(Nr, Ny) * 2 - 1);
end


function [ ] = trainReservoir(trainU, trainY)
    global W W_in W_ofb W_out noise Nr scaleU scaleY X 

    % global flags
    global extendX outputActivation reservoirStateInit
    
    dim = size(trainU,1);
    Nu  = size(trainU,2);
    assert(dim == size(trainY, 1));

    % scale trainU, trainY
    trainU = trainU ./ scaleU;
    trainY = trainY ./ scaleY;

    % initialize activations X
    
    if reservoirStateInit == 'rand'
        Xinit = tanh(10*randn(1, Nr));
    elseif reservoirStateInit == 'zero'
        Xinit = zeros(1, Nr);
    end
    X  = [Xinit; zeros(dim-1, Nr)];
    
    fprintf('#training fields: %d\n', dim);
    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end
        
    time = toc; fprintf('fitting W_out...\n')
    % compute W_out: W_out*X' = atanh(trainY)
    
    if extendX
        extX = [X, trainU];
    else
        extX = X;
    end

    % Using a pseudo inverse
    % P = pinv(X);
    % W_out = (P*atanh(trainY))';
    
    % By solving the normal equations and including Tikhonov
    % regularization
    lambda = 1e3;
    if extendX 
        Xnormal = extX'*extX + lambda * speye(Nu+Nr);
    else
        Xnormal = extX'*extX + lambda * speye(Nr);
    end
    
    if outputActivation
        f_out     = @(y) tanh(y);
        inv_f_out = @(y) atanh(y);
    else
        f_out = @(y) y;
        inv_f_out = @(y) y;
    end
    
    b       = extX'*inv_f_out(trainY);
    W_out   = (Xnormal \ b)';
     
    fprintf('fitting W_out... done (%fs)\n', toc-time)

    % get training error
    predY = f_out(extX*W_out');

    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));
end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr alpha
    pre = W*state' + W_in*u' + W_ofb*y';
    act = alpha * tanh(pre) + (1-alpha) * state' + noise * (rand(Nr,1) - 0.5);
end

function [] = plotBlockQG(x)
    global x1max x2max nxa nya
    M = numel(x);
    
    % indexing
    b1 = 1:M/2;
    b2 = (M/2)+1:M;

    % unscale
    x(b1,:) = x(b1,:) * x1max;
    x(b2,:) = x(b2,:) * x2max;
    blockid = [b1; b2];
    blockid = blockid(:);
    x = x(blockid, :);
    
    plotQG(nxa, nya, 1, x, false)
end