function [ ] = experiment(varargin)
% The core experiment is repeated with <reps>*<shifts> realisations of
% the network. The training data changes with <shifts>.

    time = tic;
    global pid procs exp_name storeState memory windowsize

    if ~isdeployed
        addpath('~/local/matlab/');
        addpath('~/Projects/ESN/matlab');
    end

    storeState = 'all';   % which states to store

    hyp = struct();
    range2str = @ (range) ['_', num2str(range(1)), '-', num2str(range(end)), '_'];

    %---------------------------------------------------------
    % settings that define the experiment
    run_pars.esn_on     = true; % enable/disable ESN
    run_pars.model_on   = true; % enable/disable equations
    run_pars.no_testing = true; % if true we ignore the test data and just run the model
    exp_id = {'ReservoirSize'};

    % numeric options
    name = 'ReservoirSize';
    hyp.(name).range   = [3000];
    hyp.(name).descr   = ['NR', range2str(hyp.(name).range)];
    hyp.(name).default = 8000;
    
    name = 'BlockSize';
    hyp.(name).range   = [1,16];
    hyp.(name).descr   = ['BS', range2str(hyp.(name).range)];
    hyp.(name).default = 8;

    name = 'TrainingSamples';
    hyp.(name).range   = [1000,2000,3000,4000,5000,6000,7000];
    hyp.(name).descr   = ['SP', range2str(hyp.(name).range)];
    hyp.(name).default = 4000;

    name = 'ReductionFactor';
    hyp.(name).range   = [16, 32];
    hyp.(name).descr   = ['RF', range2str(hyp.(name).range)];
    hyp.(name).default = 1;

    name = 'Alpha';
    hyp.(name).range   = [0.2,1.0];
    hyp.(name).descr   = ['AP', range2str(hyp.(name).range)];
    hyp.(name).default = 0.2;
    
    name = 'RhoMax';
    hyp.(name).range   = [0.3];
    hyp.(name).descr   = ['RH', range2str(hyp.(name).range)];
    hyp.(name).default = 0.3;
    
    name = 'FeedthroughAmp';
    hyp.(name).range   = [0.1,0.7,1.0];
    hyp.(name).descr   = ['FA', range2str(hyp.(name).range)];
    hyp.(name).default = 1.0;
    
    name = 'ReservoirAmp';
    hyp.(name).range   = [0.1,1.0,10];
    hyp.(name).descr   = ['RA', range2str(hyp.(name).range)];
    hyp.(name).default = 1.0;

    name = 'InAmplitude';
    hyp.(name).range   = [1.0,10.0];
    hyp.(name).descr   = ['IA', range2str(hyp.(name).range)];
    hyp.(name).default = 1.0;

    name = 'AverageDegree';
    hyp.(name).range   = [5,10,20,30];
    hyp.(name).descr   = ['AD', range2str(hyp.(name).range)];
    hyp.(name).default = 10;

    name = 'Lambda';
    hyp.(name).range   = [1e-8, 1e-6, 1e-4, 1e-1];
    hyp.(name).descr   = ['LB', range2str(hyp.(name).range)];
    hyp.(name).default = 1e-1;

    % string based options
    name = 'SquaredStates';
    hyp.(name).opts    = {'disabled', 'append', 'even'};
    hyp.(name).range   = [1, 3];
    hyp.(name).descr   = ['SS', range2str(hyp.(name).range)];
    hyp.(name).default = 3;

    name = 'ReservoirStateInit';
    hyp.(name).opts    = {'zero', 'random'};
    hyp.(name).range   = [1, 2];
    hyp.(name).descr   = ['RI', range2str(hyp.(name).range)];
    hyp.(name).default = 2;

    xlab =  exp_id;
    ylab = 'Predicted days';

    % ensemble setup
    shifts   = 1;              % shifts in training_range
    reps     = 1;              % repetitions per shift
    maxPreds = floor(50*365);  % prediction barrier (in samples)
    %---------------------------------------------------------

    % identifier -> numeric index
    hypids = fieldnames(hyp);
    id2ind = @ (str) find(strcmp(hypids, str));

    exp_ind = []; file_descr = [];
    for i = 1:numel(exp_id)
        exp_ind{i}    = id2ind(exp_id{i});
        file_descr{i} = hyp.(exp_id{i}).descr;
    end

    assert(~isempty(exp_ind));

    exp_name = [[file_descr{:}], ...
                'ESN', num2str(run_pars.esn_on), '_', ...
                'MDL', num2str(run_pars.model_on)];

    evalstr = '';

    for i = 1:numel(hypids)
        if ~sum(strcmp(hypids{i}, exp_id))
            hyp.(hypids{i}).range = hyp.(hypids{i}).default;
        end

        evalstr = [evalstr, 'hyp.(hypids{', num2str(i), '}).range'];
        if i < numel(hypids)
            evalstr = [evalstr, ', '];
        end
    end

    eval(['hyp_range = combvec(', evalstr, ');']);

    % The core experiment is repeated with <reps>*<shifts> realisations of
    % the network. The range of the training data changes with <shifts>.
    cvec = combvec((1:reps),(1:shifts))';
    rvec = cvec(:,1);
    svec = cvec(:,2);
    Ni = numel(svec); % number of indices

    num_exp       = size(hyp_range,2);
    predictions   = cell(Ni, num_exp);
    truths        = cell(Ni, num_exp);
    errs          = cell(Ni, num_exp);
    esnXsnaps     = cell(Ni, num_exp);
    num_predicted = zeros(shifts*reps, num_exp); % valid predicted time steps

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

    print0('--------------------------------------------\n')
    print0(' ----   MLQG experiment - procs  = %d \n', procs)
    print0('  ---                   - pid    = %d \n', pid)
    print0('   --   %s \n', exp_name);

    if ~exist('trdata','var')
        print0('load training data...\n'); tic;
        fname_base = 'N64-N32_ff2_Re1.0e+03-Re5.0e+02_Tstart50_Tend100';
        % fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
        trdata = load(['~/Projects/qg/matlab/MLQG/data/training/', fname_base, '.mat']);
        print0('load training data... done (%fs)\n', toc);
    end


    nxc  = trdata.nxc
    nyc  = trdata.nyc
    nun  = 2;
    dim  = nxc*nyc*nun;
    ampl = trdata.ampl
    stir = trdata.stir
    Re_c = trdata.Re_c

    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;

    % create coarse QG setup for this problem
    qgc = QG(nxc, nyc, 1);  % coarse QG with periodic bdc
    qgc.set_par(5,  Re_c);  % Reynolds number for coarse model
    qgc.set_par(11, ampl);  % stirring amplitude
    qgc.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)

    % naming change and getting stuff out of memory
    trdata.PRX = trdata.ERX;
    trdata.ERX = [];
    trdata.R   = [];
    rmfield(trdata, 'ERX');
    rmfield(trdata, 'R');

    for j = 1:num_exp
        % print experiment info
        str = [];
        for id = 1:numel(hypids)
            str{id} = ['\n    ', hypids{id}, ': ', ...
                       num2str(hyp_range(id2ind(hypids{id}), j)),...
                       ''];
        end

        str = [str, '\n'];
        print0([str{:}]);

        % set experiment parameters.
        esn_pars.Nr          = hyp_range(id2ind('ReservoirSize'), j);
        bs                   = hyp_range(id2ind('BlockSize'), j);
        samples              = hyp_range(id2ind('TrainingSamples'), j);
        RF                   = hyp_range(id2ind('ReductionFactor'), j);
        esn_pars.alpha       = hyp_range(id2ind('Alpha'), j);
        esn_pars.rhoMax      = hyp_range(id2ind('RhoMax'), j);
        esn_pars.ftAmp       = hyp_range(id2ind('FeedthroughAmp'), j);
        esn_pars.resAmp      = hyp_range(id2ind('ReservoirAmp'), j);
        esn_pars.inAmplitude = hyp_range(id2ind('InAmplitude'), j);
        esn_pars.avgDegree   = hyp_range(id2ind('AverageDegree'), j);
        esn_pars.lambda      = hyp_range(id2ind('Lambda'), j);

        esn_pars.squaredStates = ...
            hyp.SquaredStates.opts{hyp_range(id2ind('SquaredStates'), j)};
        esn_pars.reservoirStateInit = ...
            hyp.ReservoirStateInit.opts{hyp_range(id2ind('ReservoirStateInit'), j)};

        % wavelet basis
        H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
        Na = dim / RF;
        run_pars.Ha = H(1:Na,:)';
        run_pars.Hd = H(Na+1:dim,:)';
        run_pars.Na = Na;

        print0('transform input/output data with wavelet modes\n');
        trdata.HaRX  = run_pars.Ha' * trdata.RX;
        trdata.HaPRX = run_pars.Ha' * trdata.PRX;
        % print0('compute svd\n');
        % [Uwav,~,~] = svds(trdata.HaRX(:,:), 16);

        % stopping criterion returns a stopping flag based on whatever is
        % relevant for the experiment
        run_pars.stopping_criterion = @qg_stopping_criterion;

        % shifts in the training_range
        % total number of timesteps in timeseries
        Nt = size(trdata.RX, 2);
        
        % largest shift in training_range
        if run_pars.no_testing
            maxShift = Nt - samples - 1; 
        else
            maxShift = Nt - maxPreds - samples - 1; 
        end
        
        assert(maxShift > 1);
        tr_shifts = round(linspace(0, maxShift, shifts));

        % domain decomposition
        my_inds = my_indices(pid, procs, Ni);

        for i = my_inds;
            % (re)set global memory for the computation of Ke, Km, etc
            memory = struct();
            windowsize = 10;

            memory.Ha   = run_pars.Ha;
            % memory.Uwav = Uwav;

            run_pars.train_range = (1:samples)+tr_shifts(svec(i));
            run_pars.test_range  = run_pars.train_range(end) + (1:maxPreds);
            print0(' train range: %d - %d\n', ...
                    min(run_pars.train_range), max(run_pars.train_range));
            print0('  test range: %d - %d\n', ...
                    min(run_pars.test_range), max(run_pars.test_range));
            
            [predY, testY, err, esnX] = ...
                experiment_core(qgc, trdata, esn_pars, run_pars);
            

            % Run the experiment in a try/catch block, at most max_tries times.
            % try_count = 0;
            % exc_count = 0;
            % max_tries = 5;
            % while (exc_count == try_count) && (try_count < max_tries)
            %     try
            %         try_count = try_count + 1;
            %         [predY, testY, err] = ...
            %             experiment_core(qgc, trdata, esn_pars, run_pars);
            %     catch ME
            %         print0('ERROR: pid %d, i %d, %s\n', pid, i, ME.message);
            %         exc_count = exc_count + 1;
            %     end
            % end
            % if try_count >= max_tries
            %     ME = MException('experiment:fatalError', ...
            %                     'Too many fails in experiment_core...');
            %     throw(ME);
            % end

            num_predicted(i, j) = size(predY, 1);

            if strcmp(storeState, 'all')
                predictions{i, j} = predY(:,:);
                truths{i, j}      = testY(:,:);
                esnXsnaps{i,j}    = esnX(round(linspace(1,size(esnX,1),20)),:);

            elseif strcmp(storeState, 'final');
                predictions{i, j} = predY(end,:);
                truths{i, j}      = testY(end,:);
                esnXsnaps{i,j}    = esnX(end,:);
            else
                error('Unexpected input');
            end

            errs{i, j} = err;
            store_results(my_inds, hyp_range, hyp, exp_id, exp_ind, xlab, ylab, ...
                          num_predicted, errs, predictions, truths, ...
                          run_pars, esn_pars, esnXsnaps);
        end
    end
    print0('done (%fs)\n', toc(time));
end

function [] = store_results(varargin)
    global pid procs exp_name

    if procs > 1
        run_type = 'parallel';
    else
        run_type = 'serial';
    end

    path = sprintf('~/Projects/qg/matlab/MLQG/data/experiments/%s/%s', exp_name, run_type);
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
    fprintf('\n\n');
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
