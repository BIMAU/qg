%%% QG TRAINING DATA CREATION
fprintf('load full model data...\n'); tic;
%orig_base = 'N256_Re4.0e+04_Tstart151_Tend179_F0.5_Stir0_Rot1';
%orig_base = 'N128_Re1.0e+04_Tstart159_Tend187_F0.5_Stir0_Rot1'
orig_base = 'N64_Re1.0e+03_Tstart0_Tend100_F2_Stir0_Rot1';
fmdata = load(['data/fullmodel/', orig_base, '.mat']);
fprintf('load full model data... done (%fs)\n', toc);

nun = 2; % number of unknowns
nx  = fmdata.nx;
ny  = fmdata.ny;

% full model parameters
Re_f = fmdata.Re;    % Reynolds number on fine grid
ampl = fmdata.ampl;  % Forcing amplitude
stir = 0;

%% coarse model setup
ff   = 2;          % coarsing factor
Re_c = Re_f / 2;   % Reynolds number for coarse model
nxc  = nx / ff;
nyc  = ny / ff;

qgc  = QG(nxc, nyc, 1); % coarse QG with periodic bdc
qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qgc.set_par(11, ampl);  % stirring amplitude
qgc.set_par(5,  Re_c);  % Reynolds number for coarse model

%% create grid transfer matrix
% grid transfer operator for single unknown
[~,~,R,~] = gridTransfers(nx, 'periodic');
% get the correct operator for our grid ordering:
R = kron(R, [1,0;0,1]);

%% restrict states to coarse grid
Ns    = size(fmdata.states,2);
Nt    = min(Ns, 50*365); % how many days do we want?
trnge = Ns-Nt+1:Ns; 
X     = fmdata.states(:,trnge);
RX    = R*X; % restriction
fmdata.times = trnge / 365;

%% generate predictions
tpars     = {};
tpars.dt  = fmdata.dt;
tpars.Nt  = Nt;
%tpars.Nt  = 365; % a year (time step is 1 day)
tpars.T   = ceil(tpars.Nt*tpars.dt);
ERX       = generate_coarse_predictions(qgc, RX, tpars);

%% save to file
ReStr = sprintf('Re%1.1e-Re%1.1e',Re_f, Re_c);
nxStr = sprintf('N%d-N%d_ff%d', nx, nxc, ff);
tsStr = sprintf('Tstart%d', round(fmdata.times(1)));
teStr = sprintf('Tend%d',   round(fmdata.times(tpars.Nt)));

fname_base = [nxStr, '_', ReStr, '_', tsStr, '_', teStr];

fprintf('saving data to %s\n', [fname_base,'.mat']);

RX = RX(:,1:tpars.Nt);
tic
save(['data/training/', fname_base, '.mat'], ...
     'RX', 'ERX', 'orig_base', 'R', 'Re_c', ...
     'ampl', 'stir', 'nxc', 'nyc', 'tpars');
toc

function [ERX] = generate_coarse_predictions(model, RX, tpars)
    N = size(RX,1);

    % solution matrix
    ERX = zeros(N, tpars.Nt);

    time = tic;
    fprintf('Generate predictions... \n');
    avgK = 0;
    for i = 1:tpars.Nt
        if mod(i,10) == 0
            fprintf('     %d / %d \n', i, tpars.Nt);
        end
        [ERX(:,i), k] = model.step(RX(:,i), tpars.dt);
        avgK = avgK + k;
    end
    fprintf('Generate predictions... done (%f)\n', toc(time));
    fprintf('Average # Newton iterations: (%f)\n', avgK / tpars.Nt);
end