fprintf('load training data...\n'); tic;
fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
%trdata = load(['data/training/', fname_base, '.mat']);
fprintf('load training data... done (%fs)\n', toc);

nxc  = trdata.nxc;
nyc  = trdata.nyc;
nun  = 2;
dim  = nxc*nyc*nun;
RX   = trdata.RX;
ERX  = trdata.ERX;
ampl = trdata.ampl;
stir = trdata.stir;
Re_c = trdata.Re_c;

Ldim    = 1e6;
Udim    = 3.171e-2;
tdim    = Ldim / Udim; % in seconds
scaling = 3600*24/tdim;

% create coarse QG setup for this problem
qgc = QG(nxc, nyc, 1);  % coarse QG with periodic bdc
qgc.set_par(5,  Re_c);  % Reynolds number for coarse model
qgc.set_par(11, ampl);  % stirring amplitude
qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)

esn_pars.Nr = 3000;
trdata.PRX = trdata.ERX;

% dimension reduction Na
Na = dim / 2;

% create wavelet basis
bs = 32; % block size
H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
run_pars.Ha = H(1:Na,:)';
run_pars.Hd = H(Na+1:dim,:)';
run_pars.Na = Na;

run_pars.esn_on   = true;
run_pars.model_on = true;
run_pars.stopping_criterion = @qg_stopping_criterion;

cutoff = round(trdata.tpars.Nt * 3 / 4);
run_pars.train_range = 1:cutoff;
run_pars.test_range  = cutoff+1:cutoff+50;

[predY, testY] = experiment_core(qgc, trdata, esn_pars, run_pars);

% for i = 1:4:numel(run_pars.test_range)
%     plotQG(nxc,nyc,1,scaling*(predY(i,:)),false)
%     caxis([-0.25,0.25]);
%     drawnow
% end

