% Initialize full QG model
nun = 2;
nx = 256;
ny = 256;
n = nx*ny*nun;
qg = QG(nx, ny, 1);  % periodic bdc setup

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365*day;

Re   = 4e4;   % Reynolds number on fine grid
ampl = 0.5;   % Forcing amplitude
Tend = 100;   % End time, nondim. timescale is in years
dt   = 0.01;  % Time step
th   = 1.0;   % Theta
t0   = 0;     % time

% Set parameters in QG
qg.set_par(18,  0.0);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qg.set_par(11, ampl);  % stirring amplitude
qg.set_par(5,    Re);  % Reynolds number

% Also create a coarse QG model
ff  = 4;    % coarsening factor
nxa = 256 / ff;
nya = 256 / ff;
qg_a = QG(nxa, nya, 1);

% Set parameters in QG
qg_a.set_par(18, 0.0);    % stirring type: 0 = cos(5x), 1 = sin(16x)
qg_a.set_par(11, ampl);   % stirring amplitude
qg_a.set_par(5,  Re/100); % Reynolds number for coarse model

% load data: full model, fixed timestep
fprintf('load full model data...\n')
fname_base = 'N256_Re4.0e+04_Tstart141_Tend142_F0.5';
data = load(['data/fullmodel/', fname_base, '.mat']);
assert(nx == data.nx);
fprintf('load full model data... done\n')
times = data.times;

scaling = 3600*24/tdim;
crange = [-0.2,0.2];
Frange = [-10,10];

N = round(size(data.states,2));

idx = 1;
eddyF = zeros(nxa*nya*2, N);
xa    = zeros(nxa*nya*2, N);

% create averaging matrix
fprintf('building averaging operator...\n')
Pf = blockpermutation(nx,ny,nun,ff);   % full block permutation
Pa = blockpermutation(nxa,nya,nun,1);  % block permutation on
                                       % coarse grid (averages)
H  = ones(1,ff^2)/(ff^2);              %
M  = speye(n/(ff^2));                  %
H  = kron(M,H);                        % compute averages of blocks
P  = Pa'*H*Pf;                         %
fprintf('building averaging operator... done\n')

fprintf('creating eddy forcing... \n', idx, N);

tic
for idx = 1:N
    fprintf('                 %d/%d\n', idx, N);
    xf  = data.states(:,idx);    % full state
    F   = qg.rhs(xf);            % full rhs
    xa(:,idx) = P*xf;            % averaged state
    F_a = qg_a.rhs(xa(:,idx));   % coarse rhs
    eddyF(:,idx) = F_a-P*F;      % eddy forcing given by discrepancy
end
Telap = toc;
fprintf('creating eddy forcing... done: %f\n', Telap);
fname = ['eddyforcing_', fname_base, '.mat'];
fprintf('saving data to %s\n', fname);
save(['data/eddyforcing/',fname], 'eddyF', 'xa', ...
     'N', 'P', 'nxa', 'nya', 'times')