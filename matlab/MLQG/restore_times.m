% perform transients initialized with progressively less wavelet
% modes.

append  = true;
compute = true;
compute = false;

fprintf('load data...\n'); tic;
if ~exist('data','var')
    data = load('data/fullmodel/N128_Re1.0e+04_Tstart159_Tend187_F0.5_Stir0_Rot1.mat');
end
fprintf('load data... done (%fs)\n', toc);

if compute
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(4)
    end
end

% wavelet basis for full problem
bs  = 32;
nun = 2;
H   = create_wavelet_basis(data.nx, data.ny, nun, bs, true);

% initial state
n = nun*data.nx*data.ny; % total dimension
tstart = 1;
x_init = data.states(:, tstart);
assert(size(x_init,1) == n);

% orthogonal projections on reduced wavelet bases
rd = [1,1,256,512]; % reduction factors
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

totalT = size(data.states,2);
reps = 10;

T = min(200, (totalT - tstart)/reps);

ampl = [0.5, 1.5, 1.0, 1.0];
if compute
    for r = 1:reps
        if ~exist('x','var') || isempty(x) || ~append
            for j = 1:Nr
                x{j} = zeros(n,T);
                x{j}(:,1) = Pa{j}*x_init;
            end
            sidx = 2
        else
            sidx = size(x{1},2) + 1
            for j = 1:Nr
                x{j} = [x{j}, zeros(n,T)];
            end
        end

        for j = 1:Nr
            x{j} = time_integration(x{j}, dt, data, ampl(j), sidx);
        end

        fprintf('saving...\n')
        time = tic;
        savexh5(x, Nr);
        fprintf('saving... done (%f)\n', toc(time));
    end

    fprintf('computation done\n')
    return
else
    fprintf('reading...\n')
    time = tic;
    x = readxh5(Nr);
    fprintf('reading... done (%f)\n', toc(time));
end

T = size(x{1},2);

windowsize = 10;

t  = tstart+1:tstart+T;
t0 = 1:t(1)-1;

Re   = data.Re;
stir = data.stir;
ampl = data.ampl;

qg = QG(data.nx, data.ny, 1);   % QG with periodic bdc
qg.set_par(5,  Re);   % Reynolds number for coarse model
qg.set_par(18, stir); % stirring type: 0 = cos(5x), 1 = sin(16x)
qg.set_par(11, ampl); % stirring amplitude

dat = data.states(:,[t0,t]);
stats0 = get_statistics(qg, data, dat, windowsize);

dat = x{1}(:,1:T);
stats1 = get_statistics(qg, data, dat, windowsize);

dat = x{2}(:,1:T);
stats2 = get_statistics(qg, data, dat, windowsize);

dat = x{3}(:,1:T);
stats3 = get_statistics(qg, data, dat, windowsize);

dat = x{4}(:,1:T);
stats4 = get_statistics(qg, data, dat, windowsize);

% stopping criterion
% err1 = []; err2 = []; err3 = []; err4 = [];
% for k = 1:numel(t)-2;
%     test  = data.states(:,t0(end)+k);
%     pred1 = x{1}(:,k);
%     [~, err1(k), ~, ~] = qg_stopping_criterion(qg, pred1, test);
%     pred2 = x{2}(:,k);
%     [~, err2(k), ~, ~] = qg_stopping_criterion(qg, pred2, test);
%     pred3 = x{3}(:,k);
%     [~, err3(k), ~, ~] = qg_stopping_criterion(qg, pred3, test);
%     pred4 = x{4}(:,k);
%     [~, err4(k), ~, ~] = qg_stopping_criterion(qg, pred4, test);
% end

figure(1)
cols = lines(5);
subplot(3,1,1)
st = 1;
plot([t0(2:end),t(1:end-windowsize)], stats0.Km, '.-', 'markersize', 1, 'color', cols(1,:)); hold on;
plot([t0(2:end),t(1:end-windowsize)], ...
     repmat(mean(stats0.Km)+2*sqrt(var(stats0.Km)),1,numel(stats0.Km)), ...
     '--', 'markersize', 1, 'color', cols(1,:));
plot([t0(2:end),t(1:end-windowsize)], ...
     repmat(mean(stats0.Km)-2*sqrt(var(stats0.Km)),1,numel(stats0.Km)), ...
     '--', 'markersize', 1, 'color', cols(1,:));

plot(t(2:end-windowsize), stats1.Km, '.-', 'markersize', 1, 'color', cols(2,:));
plot(t(2:end-windowsize), stats2.Km, '.-', 'markersize', 1, 'color', cols(3,:));
plot(t(2:end-windowsize), stats3.Km, '.-', 'markersize', 1, 'color', cols(4,:));
plot(t(2:end-windowsize), stats4.Km, '.-', 'markersize', 1, 'color', cols(5,:)); hold off;
xlabel('days')
ylabel('Km: mean kinetic energy')
title('restore after dimension reduction')
xlim([t(2),t(end-windowsize)])

subplot(3,1,2)
f1 = plot([t0,t(2:end-windowsize)], stats0.Ke,'.-', 'markersize', 1, 'color', cols(1,:)); hold on;
plot([t0(2:end),t(1:end-windowsize)], ...
     repmat(mean(stats0.Ke)+2*sqrt(var(stats0.Ke)),1,numel(stats0.Ke)), ...
     '--', 'markersize', 1, 'color', cols(1,:));
plot([t0(2:end),t(1:end-windowsize)], ...
     repmat(mean(stats0.Ke)-2*sqrt(var(stats0.Ke)),1,numel(stats0.Ke)), ...
     '--', 'markersize', 1, 'color', cols(1,:));

f2 = plot(t(2:end-windowsize), stats1.Ke,'.-', 'markersize', 1, 'color', cols(2,:));
f3 = plot(t(2:end-windowsize), stats2.Ke,'.-', 'markersize', 1, 'color', cols(3,:));
f4 = plot(t(2:end-windowsize), stats3.Ke,'.-', 'markersize', 1, 'color', cols(4,:));
f5 = plot(t(2:end-windowsize), stats4.Ke,'.-', 'markersize', 1, 'color', cols(5,:)); hold off;
xlabel('days')
ylabel('Ke: eddy kinetic energy')
ylim([0,max(stats2.Ke)])
title(sprintf('sliding window: %d days', windowsize));
legend([f1,f2,f3,f4,f5], '1','0.5','1.5','reduced 256','reduced 512', 'location','northwest')
xlim([t(2),t(end-windowsize)])

subplot(3,1,3)
plot([t0,t(2:end-windowsize)], stats0.Zmean,'.-', 'markersize', 1, 'color', cols(1,:)); hold on;
plot([t0,t(2:end-windowsize)], ...
     repmat(mean(stats0.Zmean)+2*sqrt(var(stats0.Zmean)),1,numel(stats0.Zmean)), ...
     '--', 'markersize', 1, 'color', cols(1,:));
plot([t0,t(2:end-windowsize)], ...
     repmat(mean(stats0.Zmean)-2*sqrt(var(stats0.Zmean)),1,numel(stats0.Zmean)), ...
     '--', 'markersize', 1, 'color', cols(1,:));

plot(t(2:end-windowsize),      stats1.Zmean,'.-', 'markersize', 1, 'color', cols(2,:));
plot(t(2:end-windowsize),      stats2.Zmean,'.-', 'markersize', 1, 'color', cols(3,:));
plot(t(2:end-windowsize),      stats3.Zmean,'.-', 'markersize', 1, 'color', cols(4,:));
plot(t(2:end-windowsize),      stats4.Zmean,'.-', 'markersize', 1, 'color', cols(5,:));
hold off;
xlabel('days')
ylabel('mean enstrophy')
xlim([t(2),t(end-windowsize)])

while true
    slice = input('slice: ')

    figure(2)
    subplot(2,2,1)
    plotQG(data.nx, data.ny, 1, scaling*data.states(:, tstart+slice), false);
    caxis([-0.3,0.3]);
    subplot(2,2,2)
    plotQG(data.nx, data.ny, 1, scaling*x{1}(:,slice), false);
    caxis([-0.3,0.3]);
    subplot(2,2,3)
    plotQG(data.nx, data.ny, 1, scaling*x{2}(:,slice), false);
    caxis([-0.3,0.3]);
    subplot(2,2,4)
    plotQG(data.nx, data.ny, 1, scaling*x{4}(:,slice), false);
    caxis([-0.3,0.3]);

    figure(3)
    cols = lines(5);
    f1 = plotQGspectrum(qg, scaling*data.states(:, tstart+slice), 5, cols(1,:));
    hold on;
    f2 = plotQGspectrum(qg, scaling*x{1}(:,slice), 5, cols(2,:));
    hold on;
    f3 = plotQGspectrum(qg, scaling*x{2}(:,slice), 5, cols(3,:));
    hold on;
    f4 = plotQGspectrum(qg, scaling*x{3}(:,slice), 5, cols(5,:));
    hold on;
    f5 = plotQGspectrum(qg, scaling*x{4}(:,slice), 5, cols(4,:));
    legend([f1,f2,f3,f4,f5],'1','0.5','1.5','reduced 256','reduced 512', 'location', 'southwest')
end


%%
function [x] = time_integration(x, dt, data, ampl_fact, startidx)

    sidx = 2;
    ampl = data.ampl;
    switch nargin
      case 4
        ampl     = data.ampl * ampl_fact;
      case 5
        ampl = data.ampl * ampl_fact;
        sidx = startidx;
    end

    qg = QG(data.nx, data.ny, 1);
    qg.set_par(18, data.stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
    qg.set_par(11, ampl);  % stirring amplitude
    qg.set_par(5,  data.Re);    % Reynolds number

    % assume we start with something
    fprintf('t = %d days, ', sidx - 1 );
    fprintf(' |x| = %2.3e \n', norm(x(:,sidx-1), 2));
    startnrm = norm(x(:,sidx-1),2);
    assert(startnrm > 0);
    
    for i = sidx:size(x,2)
        fprintf('t = %d days, ', i);
        x(:,i) = qg.step(x(:,i-1),dt);
        fprintf(' |x| = %2.3e \n', norm(x(:,i), 2));
    end
end

function savexh5(x, Nr)
    [stat,cmd] = system('rm x.h5', '-echo');
    disp(cmd)
    for j = 1:Nr
        try
            h5create('x.h5', ['/',num2str(j),'/'], size(x{j}));
        catch ME
            disp(ME.message);
        end

        h5write('x.h5',  ['/',num2str(j),'/'], x{j});
    end
end

function [x] = readxh5(Nr)
    x = struct([]);
    for j = 1:Nr
        x{j} = h5read('x.h5',  ['/',num2str(j),'/']);
    end
end