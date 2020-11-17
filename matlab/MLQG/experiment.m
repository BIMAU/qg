function [ ] = experiment(pid, jobs)

    fprintf('load training data...\n'); tic;
    fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
    trdata = load(['data/training/', fname_base, '.mat']);
    fprintf('load training data... done (%fs)\n', toc);

    %# TODO parallel stuff
    pid  = 0;
    jobs = 1;

    nxc  = trdata.nxc;
    nyc  = trdata.nyc;
    nun  = 2;
    dim  = nxc*nyc*nun;
    RX   = trdata.RX;
    ERX  = trdata.ERX;
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
    trdata.PRX = trdata.ERX;

    % dimension reduction Na
    Na = dim / 2;

    % create wavelet basis
    bs = 32; % block size
    H  = create_wavelet_basis(nxc, nyc, nun, bs, true);
    run_pars.Ha = H(1:Na,:)';
    run_pars.Hd = H(Na+1:dim,:)';
    run_pars.Na = Na;

    run_pars.esn_on   = true; % enable/disable ESN
    run_pars.model_on = true; % enable/disable equations

    % stopping criterion returns a stopping flag based on whatever is
    % relevant for the experiment
    run_pars.stopping_criterion = @qg_stopping_criterion;

    samples   = 3000;    % samples in the training_range
    shifts    = 5;       % shifts in training_range
    maxShift  = 4000;    % largest shift in training_range
    reps      = 10;      % repetitions
    preds     = zeros(shifts*reps,1); % predictions
    maxPreds  = 365;
    tr_shifts = round(linspace(0,maxShift,shifts)); % shifts in the training_range

    cvec = combvec((1:reps),(1:shifts))';
    rvec = cvec(:,1);
    svec = cvec(:,2);

    Ni = numel(svec);
    my_inds = my_indices(pid, jobs, Ni);

    for i = my_inds;
        train_range = (1:samples)+tr_shifts(svec(i));
        fprintf(' train range: %d - %d\n', min(train_range), max(train_range));
        run_pars.train_range = train_range;
        run_pars.test_range  = train_range(end) + (1:maxPreds);
        [predY, testY] = experiment_core(qgc, trdata, esn_pars, run_pars);
        preds(i) = size(predY, 1);
    end

end

function [inds] = my_indices(pid, jobs, Ni)
% a simple decomposition to take care of parallel needs

    assert((pid < jobs) && (pid >= 0))

    k  = jobs
    decomp = [];
    offset = 0;
    remain = Ni; % elements that remain
    decomp = cell(jobs);
    for i = 1:jobs
        subset    = floor(remain / k);
        decomp{i} = offset + 1: offset + subset;
        offset    = offset + subset;

        k      = k-1;
        remain = remain - subset;
    end

    inds = decomp{pid+1,:};
end
