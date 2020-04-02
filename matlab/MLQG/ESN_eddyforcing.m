global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY Rstate

% load eddy forcing and coarse state data
fname = 'eddyforcing_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
data  = load(['data/eddyforcing/', fname, '.mat']);

% setup input data U: averages and eddy forcings
U  = [data.xa(:,1:end-1); data.eddyF(:,1:end-1)]';
Nu = size(U,2);

% setup output data Y: eddy forcings, one (fixed) timestep further
Y  = [data.eddyF(:,2:end)'];
Ny = size(Y,2);

createReservoir();

trainReservoir(U, Y);

for i = 2:size(Y,1)

    subplot(1,2,1)
    plotQG(nxa,nya,1,scaling*data.states(:,idx),false)
    caxis(crange);
    titleString = sprintf('Fine vorticity (day^{-1}), t = %3.0fd', ...
                          (data.times(idx)-data.times(1)) / day);
    title(titleString);

    subplot(1,2,2)
    plotQG(nxa, nya, 1, scaling*eddyF(:,idx), false);
    colorbar
    caxis(Frange);
    title('Eddy forcing for coarse model')
    drawnow

end


function [ ] = createReservoir()
    global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY Rstate

    % Reservoir parameters
    Nr       = 300;
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
    scaleU = 0.99*max(abs(trainU(:)));
    trainU = trainU / scaleU;
    scaleY = 0.99*max(abs(trainY(:)));
    trainY = trainY / scaleY;

    % initialize activations X
    X = zeros(dim, Nr);

    % iterate the state, save all neuron activations in X
    for k = 2:dim
        X(k, :) = update(X(k-1, :), trainU(k, :), trainY(k-1, :));
    end

    % Extend X with raw input columns
    %extX = [X, trainU];
    extX = X;

    % Create pseudo inverse
    P = pinv(extX);

    % compute W_out
    W_out = (P*atanh(trainY))';
    predY = tanh(extX*W_out');
    fprintf('#training fields: %d\n', dim);
    fprintf('training error: %e\n', sqrt(mean((predY(:) - trainY(:)).^2)));

    % save last reservoir state to use in rhs
    Rstate = X(end,:);
end

function [act] = update(state, u, y)
    global W W_in W_ofb noise Nr
    pre = W*state' + W_in*u' + W_ofb*y';
    act = tanh(pre) + noise * (rand(Nr,1) - 0.5);
end