fprintf('load training data...\n'); tic;
%fname_base = 'N256-N128_ff2_Re4.0e+04-Re4.0e+02_Tstart151_Tend152';
%trdata = load(['data/training/', fname_base, '.mat']);
fprintf('load training data... done (%fs)\n', toc);

nxc = trdata.nxc;
nyc = trdata.nyc;
nun = 2;
dim = nxc*nyc*nun;
RX  = trdata.RX;
ERX = trdata.ERX;
%ampl = trdata.ampl;
ampl = 0.5;
%stir = trdata.stir;
stir = 0;
Re_c = trdata.Re_c;

% create coarse QG setup for this problem
qgc = QG(nxc, nyc, 1);  % coarse QG with periodic bdc
qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qgc.set_par(11, ampl);  % stirring amplitude
qgc.set_par(5,  Re_c);  % Reynolds number for coarse model

% create wavelet basis, set block size
bs = 32;
H  = create_wavelet_basis(nxc, nyc, nun, bs, true);

% dimension reduction Na
Na = dim / 16;
Ha = H(1:Na,:)';
Hd = H(Na+1:dim,:)';

% create input/output data for hybrid ESN
U = [Ha' * RX(:,1:end-1); Ha' * ERX(:,1:end-1)];
Y = [Ha' * RX(:,2:end )];

% separate training and testing data
cutoff = round(trdata.tpars.Nt * 3 / 4)
fprintf(' training time: %d\n', cutoff)
fprintf(' test time: %d\n', trdata.tpars.Nt - cutoff)
trainU = U(:,1:cutoff)';
trainY = Y(:,1:cutoff)';
testU  = U(:,cutoff+1:end)';
fullY  = RX(:,2:end);
fullY  = fullY(:,cutoff+1:end);
testY  = fullY';
% testY  = Y(:,cutoff+1:end)';

% ESN parameters
pars                    = {};
pars.scalingType        = 'standardize';
pars.Nr                 = 3000;
pars.rhoMax             = 0.3;
pars.alpha              = 1.0;
pars.Wconstruction      = 'avgDegree';
pars.avgDegree          = 30;
pars.lambda             = 1e-1;
pars.bias               = 0.0;
pars.squaredStates      = 'even';
pars.reservoirStateInit = 'random';
pars.inputMatrixType    = 'balancedSparse';
pars.inAmplitude        = 1.0;
pars.feedThrough        = true;
pars.ftRange            = Na+1:2*Na;

esn = ESN(pars.Nr, size(U,1), size(Y,1));
esn.setPars(pars);
esn.initialize;
esn.train(trainU, trainY);

Npred = size(testU, 1);
state = esn.X(end,:);
predY = zeros(Npred, dim);
predS = zeros(Npred, dim);

yk  = RX(:,cutoff+1);
yks = RX(:,cutoff+1);

for i = 1:Npred
    fprintf('hybrid reservoir step %d\n', i);
    % coarse model prediction:
    Eyk   = qgc.step(yk, trdata.tpars.dt);
    
    % append prediction to regular input data
    u_in  = [Ha'*yk(:); Ha'*Eyk(:)]';
    u_in  = esn.scaleInput(u_in);
    state = esn.update(state, u_in)';
    u_out = esn.apply(state, u_in);
    u_out = esn.unscaleOutput(u_out);
    yk    = Ha*u_out(:) + Hd*(Hd'*Eyk);

    predY(i,:) = yk;
    
    % standalone step --- 
    fprintf('      standalone step %d\n', i);
    yks = qgc.step(yks, trdata.tpars.dt);
    predS(i,:) = yks;
    
    
    
    subplot(2,2,1)
    plotQG(nxc,nyc,1,testY(i,:),false)

    subplot(2,2,2)
    vecnrmY  = vecnorm(testY(1:i,:),2,2);
    diff_hyb = (testY(1:i,:) - predY(1:i,:)) ./ vecnrmY;
    diff_std = (testY(1:i,:) - predS(1:i,:)) ./ vecnrmY;

    plot(vecnorm(diff_hyb,2,2)); hold on
    plot(vecnorm(diff_std,2,2)); hold off

    subplot(2,2,3)
    plotQG(nxc,nyc,1,yks,false)

    subplot(2,2,4)
    plotQG(nxc,nyc,1,yk,false)

    drawnow

end