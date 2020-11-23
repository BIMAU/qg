function [ ] = experiment(varargin)
% The core experiment is repeated with <reps>*<shifts> realisations of
% the network. The training data changes with <shifts>.
    time = tic;
    
    global pid procs exp_name storeState

    if ~isdeployed
        addpath('~/local/matlab/');
        addpath('~/Projects/ESN/matlab');
    end

    switch nargin
      case 0
        pid   = 0;
        procs = 1;
      case 2
        pid   = str2num(varargin{1}); % assuming input is a string
        procs = str2num(varargin{2});
      otherwise
        error('Unexpected input');
    end

    % parallel seed
    tm = clock;
    rng(round(100*pid*sqrt(tm(end))));

    exp_name   = 'test';   % experiment name
    storeState = 'final';  % which states to store
    
    fprintf('--------------------------------------------\n')
    fprintf(' ----   MLQG experiment - procs  = %d \n', procs)
    fprintf('  ---                   - pid    = %d \n', pid)
    fprintf('   --   %s \n', exp_name);

    fprintf('load training data...\n'); tic;
    fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
    trdata = load(['data/training/', fname_base, '.mat']);
    fprintf('load training data... done (%fs)\n', toc);

    nxc  = trdata.nxc;
    nyc  = trdata.nyc;
    nun  = 2;
    dim  = nxc*nyc*nun;
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

    % naming change
    trdata.PRX = trdata.ERX;
    rmfield(trdata, 'ERX');

    % dimension reduction Na
    Na = dim / 2;

    % create wavelet basis
    bs = 32; % block size
    H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
    run_pars.Ha = H(1:Na,:)';
    run_pars.Hd = H(Na+1:dim,:)';
    run_pars.Na = Na;

    fprintf('transform input/output data with wavelet modes\n');
    trdata.HaRX  = run_pars.Ha' * trdata.RX;
    trdata.HaPRX = run_pars.Ha' * trdata.PRX;
    rmfield(trdata, 'PRX');  %we do not need this one

    run_pars.esn_on   = true; % enable/disable ESN
    run_pars.model_on = true; % enable/disable equations

    % stopping criterion returns a stopping flag based on whatever is
    % relevant for the experiment
    run_pars.stopping_criterion = @qg_stopping_criterion;

    samples    = 3000;    % samples in the training_range
    shifts     = 12;      % shifts in training_range
    maxShift   = 4000;    % largest shift in training_range
    reps       = 8;       % repetitions
    maxPreds   = 365;
    tr_shifts  = round(linspace(0, maxShift, shifts)); % shifts in the training_range

    % The core experiment is repeated with <reps>*<shifts> realisations of
    % the network. The range of the training data changes with <shifts>.

    cvec = combvec((1:reps),(1:shifts))';
    rvec = cvec(:,1);
    svec = cvec(:,2);

    Ni = numel(svec);
    my_inds = my_indices(pid, procs, Ni);

    predictions   = cell(Ni,1);
    truths        = cell(Ni,1);
    errs          = cell(Ni,1);
    num_predicted = zeros(shifts*reps, 1); % valid predicted time steps

    for i = my_inds;
        train_range = (1:samples)+tr_shifts(svec(i));
        fprintf(' train range: %d - %d\n', min(train_range), max(train_range));

        run_pars.train_range = train_range;
        run_pars.test_range  = train_range(end) + (1:maxPreds);
        [predY, testY, err] = ...
            experiment_core(qgc, trdata, esn_pars, run_pars);

        num_predicted(i) = size(predY, 1);

        if strcmp(storeState, 'all')
            predictions{i} = predY(:,:);
            truths{i} = testY(:,:);
        elseif strcmp(storeState, 'final');
            predictions{i} = predY(end,:);
            truths{i} = testY(end,:);
        end

        errs{i} = err;
        store_results(my_inds, num_predicted, errs, predictions, truths);
    end
    fprintf('done (%fs)\n', toc(time));
end

function [] = store_results(my_inds, num_predicted, errs, predictions, truths)
    global pid procs exp_name

    if procs > 1
        run_type = 'parallel';
    else
        run_type = 'serial';
    end

    path = sprintf('data/experiments/%s/%s', exp_name, run_type);
    syscall = sprintf('mkdir -p %s', path);
    system(syscall);

    if strcmp(run_type, 'parallel')
        fname = sprintf('%s/results_%d.mat', path, pid);
    elseif strcmp(run_type, 'serial')
        fname = sprintf('%s/results.mat', path);
    end

    fprintf('saving results to %s\n', fname);
    save(fname, 'my_inds', 'num_predicted', 'errs', 'predictions', 'truths');
end

function [inds] = my_indices(pid, procs, Ni)
% a simple decomposition to take care of parallel needs

    assert((pid < procs) && (pid >= 0), ...
           ['assertion failed, pid ', num2str(pid)])

    k = procs;
    decomp = [];
    offset = 0;
    remain = Ni; % elements that remain
    decomp = cell(procs);
    for i = 1:procs
        subset    = floor(remain / k);
        decomp{i} = offset + 1: offset + subset;
        offset    = offset + subset;

        k      = k-1;
        remain = remain - subset;
    end

    inds = decomp{pid+1,:};
end