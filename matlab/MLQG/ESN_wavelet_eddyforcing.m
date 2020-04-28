rng(7)
global W W_in W_ofb W_out noise Nr Nu Ny scaleU scaleY X

% QG parameters:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;


% Load eddy forcing and coarse state data
% fname = 'eddyforcing_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart141_Tend142_F0.5';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend169_F0.5_Stir0_Rot1';
% fname = 'eddyforcing_ff4_Re4.0e+02_N256_Re4.0e+04_Tstart142_Tend151_F0.5_Stir0_Rot1';
% data  = load(['data/eddyforcing/', fname, '.mat']);
Nu    = size(data.xa,1);
assert(Nu == size(data.eddyF,1));

nxa = data.nxa;
nya = data.nya;
nun = 2;
assert(Nu == nxa*nya*nun);

bs = 4; % blocksize
P  = blockpermutation(nxa,nya,nun,bs);

H = haarmat(bs^2);
M = speye(Nu/(bs^2));
H = kron(M,H);
H = (H*P)';

% so now we have u = H_A*u_A + H_D*u_D with H = H_A truncated

r    = data.eddyF(:,1);
rbar = H'*r;
[~,id] = sort(abs(rbar));
dim = 200;
id = id(end-dim+1:end);

%[~,id] = sort(abs(ubar));

%id = id(end-dim+1:end);
H = H(:,id);     % truncate modes and coefficients
rbar = rbar(id); 

test = zeros(dim,1);
test(10) = 1;
r = H*rbar;
subplot(2,2,2);
plotQG(nxa, nya, 1, day*r, false)
colorbar