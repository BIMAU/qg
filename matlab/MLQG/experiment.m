function [ ] = experiment(varargin)
% The core experiment is repeated with <reps>*<shifts> realisations of
% the network. The training data changes with <shifts>.
    time = tic;
    global pid procs exp_name storeState

    if ~isdeployed
        addpath('~/local/matlab/');
        addpath('~/Projects/ESN/matlab');
    end
    
    exp_name   = 'fulldimNr10000-12000';  % experiment name
    storeState = 'final';    % which states to store

    % hyperparameter range
    hyp_range  = [10000, 12000];
    xlab       = 'Nr';
    ylab       = 'Predicted days';
 
    switch nargin
      case 0
        pid   = 0;
        procs = 1;

      case 2
        pid   = arg2value(varargin{1});
        procs = arg2value(varargin{2});
        
      case 3 
        pid    = arg2value(varargin{1});
        procs  = arg2value(varargin{2});
        trdata = varargin{3};
        
      otherwise
        error('Unexpected input');
    end

    % parallel seed
    tm = clock;
    rng(round(100*pid*sqrt(tm(end))));

    fprintf('--------------------------------------------\n')
    fprintf(' ----   MLQG experiment - procs  = %d \n', procs)
    fprintf('  ---                   - pid    = %d \n', pid)
    fprintf('   --   %s \n', exp_name);

    if ~exist('trdata','var')
        fprintf('load training data...\n'); tic;
        fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
        trdata = load(['data/training/', fname_base, '.mat']);
        fprintf('load training data... done (%fs)\n', toc);
    end

    nxc  = trdata.nxc;
    nyc  = trdata.nyc;
    nun  = 2;
    dim  = nxc*nyc*nun;
    ampl = trdata.ampl;
    stir = trdata.stir;
    Re_c = trdata.Re_c;

    % dimension reduction Na
    Na = dim;

    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;

    % create coarse QG setup for this problem
    qgc = QG(nxc, nyc, 1);  % coarse QG with periodic bdc
    qgc.set_par(5,  Re_c);  % Reynolds number for coarse model
    qgc.set_par(11, ampl);  % stirring amplitude
    qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)

    % naming change
    trdata.PRX = trdata.ERX;
    rmfield(trdata, 'ERX'); % we do not need this field, save some memory

    % create wavelet basis
    bs = 32; % block size
    H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
    run_pars.Ha = H(1:Na,:)';
    run_pars.Hd = H(Na+1:dim,:)';
    run_pars.Na = Na;

    fprintf('transform input/output data with wavelet modes\n');
    trdata.HaRX  = run_pars.Ha' * trdata.RX;
    trdata.HaPRX = run_pars.Ha' * trdata.PRX;
    rmfield(trdata, 'PRX');  % we do not need this field, save some memory

    run_pars.esn_on   = true; % enable/disable ESN
    run_pars.model_on = true; % enable/disable equations

    % stopping criterion returns a stopping flag based on whatever is
    % relevant for the experiment
    run_pars.stopping_criterion = @qg_stopping_criterion;

    samples    = 3000;    % samples in the training_range
    shifts     = 25;      % shifts in training_range
    maxShift   = 4000;    % largest shift in training_range
    reps       = 4;       % repetitions
    maxPreds   = 365;
    tr_shifts  = round(linspace(0, maxShift, shifts)); % shifts in the training_range

    num_trials = numel(hyp_range);

    % The core experiment is repeated with <reps>*<shifts> realisations of
    % the network. The range of the training data changes with <shifts>.
    cvec = combvec((1:reps),(1:shifts))';
    rvec = cvec(:,1);
    svec = cvec(:,2);

    Ni = numel(svec);
    my_inds = my_indices(pid, procs, Ni);

    predictions   = cell(Ni, num_trials);
    truths        = cell(Ni, num_trials);
    errs          = cell(Ni, num_trials);
    num_predicted = zeros(shifts*reps, num_trials); % valid predicted time steps

    for j = 1:num_trials
        esn_pars.Nr = hyp_range(j);
        
        for i = my_inds;
            run_pars.train_range = (1:samples)+tr_shifts(svec(i));
            run_pars.test_range  = run_pars.train_range(end) + (1:maxPreds);
            fprintf(' train range: %d - %d\n', min(run_pars.train_range), max(run_pars.train_range));
            fprintf('  test range: %d - %d\n', min(run_pars.test_range), max(run_pars.test_range));
            
            % Run the experiment in a try/catch block, at most max_tries times.
            % Sometimes it runs out of memory.
            try_count = 0;
            exc_count = 0;
            max_tries = 5;
            while (exc_count == try_count) && (try_count < max_tries)
                try
                    try_count = try_count + 1;
                    [predY, testY, err] = ...
                        experiment_core(qgc, trdata, esn_pars, run_pars);
                catch ME
                    fprintf('ERROR: pid %d, i %d, %s\n', pid, i, ME.message);
                    exc_count = exc_count + 1;
                end                
            end
            if try_count >= max_tries
                ME = MException('experiment:fatalError', 'Too many fails in experiment_core...');
                throw(ME);
            end
                
            num_predicted(i, j) = size(predY, 1);
            
            if strcmp(storeState, 'all')
                predictions{i, j} = predY(:,:);
                truths{i, j} = testY(:,:);

            elseif strcmp(storeState, 'final');
                predictions{i, j} = predY(end,:);
                truths{i, j} = testY(end,:);
            end

            errs{i, j} = err;
            store_results(my_inds, hyp_range, xlab, ylab, ...
                          num_predicted, errs, predictions, truths);
        end
    end
    fprintf('done (%fs)\n', toc(time));
end

function [] = store_results(varargin)
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
    for i = 1:nargin
        var = inputname(i);
        eval([var, '= varargin{i};']); 
        fprintf(' %s', var);
        if (i == 1)
            save(fname, var);
        else
            save(fname, var, '-append');
        end            
    end
    fprintf('\n');
end

function [inds] = my_indices(pid, procs, Ni)
% a simple decomposition to take care of parallel needs

    assert((pid < procs) && (pid >= 0), ...
           ['assertion erred, pid ', num2str(pid)])

    assert((procs <= Ni), ...
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