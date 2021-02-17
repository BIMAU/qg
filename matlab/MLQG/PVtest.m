fprintf('load data...\n'); tic;
if ~exist('data','var')
    data = load('data/training/N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187.mat');
    % data = load('data/fullmodel/N128_Re1.0e+04_Tstart159_Tend187_F0.5_Stir0_Rot1.mat');
end
fprintf('load data... done (%fs)\n', toc);

nx = data.nxc;
ny = data.nyc;
am = data.ampl;
Re = data.Re_c;
st = data.stir;
dt = data.tpars.dt;

nun  = 2;
m    = ny;
n    = nx;
ndim = nun*m*n;
adim = 25*ndim;
hdim = 6e+02;
f0dim = 1.0e-04;
beta0dim = 1.6e-11;
Lxdim = 1.0e+06;
Lydim = 1.0e+06;
udim = 3.171e-02;    % 3.171e-02 gives timescale of ~1year.
                     % udim = 1.6e-2;
gdim = 2.0e-02;
rhodim = 1.0e+03;
taudim = 1.5e-01;
Ahdim = 1.0e+03;
bfdim = 0.0e+00;
xmin = 0.0;
xmax = 1.0;
ymin = 0.0;
ymax = Lydim/Lxdim;

beta = beta0dim*Lxdim*Lxdim/udim

Ldim       = 1e6;
Udim       = 3.171e-2;
tdim       = Ldim / Udim;    % in seconds
om_scaling = 3600*24 / tdim; % per day

qg = QG(nx, ny, 1);  % coarse QG with periodic bdc
qg.set_par(5,  Re);  % Reynolds number for coarse model
qg.set_par(11, am);  % stirring amplitude
qg.set_par(18, st);  % stirring type: 0 = cos(5x), 1 = sin(16x)

slice = 6184;
x  = data.RX(:,slice);
x0 = data.RX(:,slice-1);

om  = reshape(x(1:2:end),  nx, ny);
ps  = reshape(x(2:2:end),  nx, ny);
om0 = reshape(x0(1:2:end), nx, ny);
ps0 = reshape(x0(2:2:end), nx, ny);

[u, v]     = qg.compute_uv(x(:));
[omx, omy] = qg.gradient(om(:));

subplot(2,2,1)
imagesc(reshape(om, nx, ny)'); colorbar

subplot(2,2,2)
imagesc(reshape(ps, nx, ny)'); colorbar

subplot(2,2,3)
ddtom = (om(:)-om0(:))/dt;
imagesc(reshape(ddtom + u.*omx+v.*omy, nx, ny)'); colorbar

ddtPV = ddtom + u.*omx+v.*omy+beta*v;

subplot(2,2,4)
imagesc(reshape(ddtPV, nx, ny)'); colorbar

mass  = [];
ddtPV = [];
for i = 2:10000
    x   = data.RX(:,i);
    x0  = data.RX(:,i);
    om  = x(1:2:end);
    om0 = x0(1:2:end);
    ps  = x(2:2:end);
    ps0 = x0(2:2:end);
    
    mass(i) = sum(ps(:));
    
    [u, v]     = qg.compute_uv(x(:));
    [omx, omy] = qg.gradient(om(:));
    ddtom = (om(:)-om0(:))/dt;
    ddtPV(i) = sum(ddtom + u.*omx+v.*omy+beta*v);
end

subplot(2,2,1)
plot(ddtPV);