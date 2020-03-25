% Initialize QG
nun = 2;
nx = 256;
ny = 256;
n = nx*ny*nun;
qg = QG(nx, ny, 1);

% QG params:
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

% Also create a coarse variant
ff  = 2; % coarsening factor
nxc = 256 / ff;
nyc = 256 / ff;
qg_a = QG(nxc, nyc, 1);

% Set parameters in QG
qg_a.set_par(18,   0.0);    % stirring type: 0 = cos(5x), 1 = sin(16x)
qg_a.set_par(11,   ampl);   % stirring amplitude
qg_a.set_par(5,    Re/100);  % Reynolds number

% load data for the fine setup 
% data = load('N256_Re4.0e+04_Tstart141_Tend142_F0.5.mat');
assert(nx == data.nx);

scaling = 3600*24/tdim;
crange = [-0.2,0.2];
Frange = [-10,10];

N = size(data.states,2);

idx = 1;

eddyF = zeros(nxc*nyc*2, N);

% create averaging matrix
fprintf('building averaging operator...\n')
Pf = blockpermutation(nx,ny,nun,ff);   % full block permutation
Pc = blockpermutation(nxc,nyc,nun,1);  % coarse block permutation
H  = ones(1,ff^2)/(ff^2);
M  = speye(n/(ff^2));
H  = kron(M,H);                        % compute averages of blocks
P  = Pc'*H*Pf;                         
fprintf('building averaging operator... done\n')

tic
for idx = 1:10
    %fprintf('creating eddy forcing... %d/%d\n', idx, N);           
    %xf  = data.states(:,idx);    
    %F = qg.rhs(xf);
    %F = average(F,nx,ny,nun,ff); 
    %F_a = qg_a.rhs(average(xf,nx,ny,nun,ff));
    %eddyF(:,idx) = F_a-F;
    
    fprintf('creating eddy forcing... %d/%d\n', idx, N);           
    xf  = data.states(:,idx);    
    F = qg.rhs(xf);
    F_a = qg_a.rhs(P*xf);
    eddyF(:,idx) = F_a-P*F;

    subplot(1,2,1)
    plotQG(nx,ny,1,scaling*xf,false)
    caxis(crange);
    titleString = sprintf('Fine vorticity (day^{-1}), t = %3.0fd', ...
                          (data.times(idx)-data.times(1)) / day);
    title(titleString);
    
    subplot(1,2,2)
    plotQG(nxc, nyc, 1, scaling*eddyF(:,idx), false);
    colorbar
    caxis(Frange);
    title('Eddy forcing for coarse model')
    
end
toc