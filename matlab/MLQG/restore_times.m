time = tic;
% perform transients initialized with progressively less wavelet
% modes.

fprintf('load data...\n'); tic;
if ~exist('data','var')
    data = load('data/fullmodel/N128_Re1.0e+04_Tstart159_Tend187_F0.5_Stir0_Rot1.mat');
end
fprintf('load data... done (%fs)\n', toc);

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(3)
end

% wavelet basis for full problem
bs  = 32;
nun = 2;
H   = create_wavelet_basis(data.nx, data.ny, nun, bs, true);

% initial state
n  = nun*data.nx*data.ny; % total dimension
x_init = data.states(:,end);
assert(size(x_init,1)==n);

% orthogonal projections on reduced wavelet bases
rd = [1,4,8,16,32,64]; % reduction factors
Nr = numel(rd);
Pa = cell(Nr,1);
Pd = cell(Nr,1);
for k = 1:Nr;
    Na    = n / rd(k);
    Ha    = H(1:Na,:)';
    Hd    = H(Na+1:n,:)';
    Pa{k} = Ha*Ha';
    Pd{k} = Hd*Hd';
end
Ldim    = 1e6;
Udim    = 3.171e-2;
tdim    = Ldim / Udim; % in seconds
scaling = 3600*24/tdim;

dt = data.dt;

T = 100;

for k = 1:Nr
    x{k} = zeros(n,T);
    x{k}(:,1) = Pa{k}*x_init - Pd{k}*x_init;
end

parfor j = 1:Nr
    x{j} = time_integration(x{j}, dt, data)
end

toc(time)
return
%%
cols = lines(10);
S1 = zeros(T,1);
for i = 1:1:T
    subplot(2,3,4)
    [f1,S1(i)]   = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{1}(:,i),5,cols(1,:));
    hold on
    [f4,S4(i)]   = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{2}(:,i),5,cols(2,:));
    hold on
    [f8,S8(i)]   = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{3}(:,i),5,cols(3,:));
    hold on
    [f16,S16(i)] = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{4}(:,i),5,cols(4,:));
    hold on
    [f32,S32(i)] = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{5}(:,i),5,cols(5,:));
    hold on
    [f64,S64(i)] = plotQGspectrum(qg_orig, data.nx, data.ny, scaling*x{6}(:,i),5,cols(6,:));
    hold off     

    legend([f1,f4,f8,f16,f32,f64], 'full', '4', '8', '16', '32', '64','location','southwest');
    xlim([20,90])
    ylim([0,0.78])

    subplot(2,3,1)
    plotQG(data.nx,data.ny,1,scaling*x{1}(:,i),false);
    title(sprintf('Na = n, t = %d days\n', i));

    subplot(2,3,2)
    plotQG(data.nx,data.ny,1,scaling*x{4}(:,i),false);
    title(sprintf('Na = n/16, t = %d days\n', i));

    subplot(2,3,3)
    plotQG(data.nx,data.ny,1,scaling*x{6}(:,i),false);
    title(sprintf('Na = n/64, t = %d days\n', i));

    drawnow
    input('')
end



%%
figure
plot(S1(2:end));
hold on
plot(S4(2:end));
plot(S8(2:end));
plot(S16(2:end));
plot(S32(2:end));
plot(S64(2:end));

hold off


%%
function [x] = time_integration(x, dt, data)
    qg = QG(data.nx, data.ny, 1);
    qg.set_par(18, data.stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
    qg.set_par(11, data.ampl);  % stirring amplitude
    qg.set_par(5,  data.Re);    % Reynolds number

    for i = 2:size(x,2)
        fprintf('t = %d days\n', i);
        x(:,i) = qg.step(x(:,i-1),dt);
    end
end