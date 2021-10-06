function [predY, testY, err, esnX] = experiment_core(model, tr_data, esn_pars, run_pars)
    % Core routine for running an experiment
    % model:         qg or ks
    % tr_data:       consists of a set of restricted 'ground truth' data
    %                and coarse model predictions based on restricted
    %                ground truth data.
    % predY:         full dimensional predictions
    % testY:         full dimensional truths

    dt  = tr_data.tpars.dt;   % time step size
    dim = size(tr_data.RX,1); % sample/state vector dimension

    % additional dimension reduction is specified in run_pars
    Na = run_pars.Na;
    Ha = run_pars.Ha;
    Hd = run_pars.Hd;

    % three different experiments based on settings in run_pars
    hybrid     = logical( (run_pars.esn_on   == true ) && ...
                          (run_pars.model_on == true ) );
    esn_only   = logical( (run_pars.esn_on   == true ) && ...
                          (run_pars.model_on == false) );
    model_only = logical( (run_pars.esn_on   == false) && ...
                          (run_pars.model_on == true ) );

    % create input/output data for hybrid or standalone ESN
    if hybrid
        print0('Create input/output data for hybrid ESN\n');
        exp_type = 'hybrid';
        if isfield(tr_data, 'HaRX') && isfield(tr_data, 'HaPRX')
            U = [tr_data.HaRX(:,1:end-1); tr_data.HaPRX(:,1:end-1)];
            Y = [tr_data.HaRX(:,2:end )];
        else
            U = [Ha' * tr_data.RX(:, 1:end-1); Ha' * tr_data.PRX(:, 1:end-1)];
            Y = [Ha' * tr_data.RX(:, 2:end )];
        end
    elseif esn_only
        print0('Create input/output data for standalone ESN\n');
        exp_type = 'esn_only';
        if isfield(tr_data, 'HaRX') && isfield(tr_data, 'HaPRX')
            U = [tr_data.HaRX(:,1:end-1)];
            Y = [tr_data.HaRX(:,2:end)];
        else
            U = [Ha' * tr_data.RX(:,1:end-1)];
            Y = [Ha' * tr_data.RX(:,2:end)];
        end
    elseif model_only
        exp_type = 'model_only';
    end

    if run_pars.esn_on
        % set training and testing data
        assert(run_pars.train_range(end)+1 == run_pars.test_range(1));
        trainU = U(:, run_pars.train_range)'; % input training
        trainY = Y(:, run_pars.train_range)'; % output training
        
        % clean up
        clear U Y
    end
    
    % full dimensional output testing
    if run_pars.no_testing
        testY = [];
    else
        testY = tr_data.RX(:, 2:end);
        testY = testY(:, run_pars.test_range)';
    end

    Npred = numel(run_pars.test_range); % number of prediction steps
    predY = zeros(Npred, dim); % full dimensional predictions
    err   = zeros(Npred, 1);   % error array
    esnX  = [];

    % initialization for the predictions
    init_idx = run_pars.train_range(end)+1;
    yk = tr_data.RX(:, init_idx);
    
    print0('initialization index: %d\n', init_idx);
    
    % clean up
    clear tr_data

    if run_pars.esn_on
        % set ESN parameters
        esn_pars = default_esn_parameters(esn_pars);

        % enable hybrid input design
        if hybrid
            esn_pars.feedThrough = true;
            esn_pars.ftRange     = Na+1:2*Na;
        end

        % create ESN, train the ESN and save the final state
        esn = ESN(esn_pars.Nr, size(trainU,2), size(trainY,2));
        esn.setPars(esn_pars);
        esn.initialize;
        esn.train(trainU, trainY);
        esn_state = esn.X(end,:);
        esnX = esn.X;
    end

    clear trainU trainY

    verbosity = 100;
    for i = 1:Npred
        
        % model prediction of next time step
        Pyk = model.step(yk, dt);

        if model_only
            % result is not adapted
            yk = Pyk;
        else
            % create an input vector for the ESN
            if hybrid
                u_in  = [Ha' * yk(:); Ha' * Pyk(:)]';
            elseif esn_only
                u_in  = [Ha' * yk(:)]';
            else
                error('incorrect model config');
            end

            u_in      = esn.scaleInput(u_in);
            esn_state = esn.update(esn_state, u_in)';
            u_out     = esn.apply(esn_state, u_in);
            u_out     = esn.unscaleOutput(u_out);

            % combine ESN prediction with model prediction
            yk = Ha*u_out(:) + Hd*(Hd'*Pyk);
        end

        predY(i,:) = yk;

        % check stopping criterion
        stop = false;
        if ~run_pars.no_testing
            [stop, err(i)] = run_pars.stopping_criterion(model, predY(i,:), testY(i,:));
        end

        if stop
            break;
        end
        
        if mod(i,verbosity) == 0
            print0('     prediction step %4d/%4d, error %1.2e, %s\n', ...
                   i, Npred, err(i), exp_type);
        end
    end
    print0(' ~~~ last prediction step %4d/%4d, error %1.2e, %s\n', ...
           i, Npred, err(i), exp_type);
    
    % truncate output arrays
    predY = predY(1:i,:);
    if ~run_pars.no_testing
        testY = testY(1:i,:);
    end
    err = err(1:i);    
end

function [pars_out] = default_esn_parameters(pars_in)
    pars_out                    = {};
    pars_out.scalingType        = 'standardize';
    pars_out.Nr                 = 3000;
    pars_out.rhoMax             = 0.3;
    pars_out.alpha              = 1.0;
    pars_out.Wconstruction      = 'avgDegree';
    pars_out.avgDegree          = 10;
    pars_out.lambda             = 1e-1;
    pars_out.bias               = 0.0;
    pars_out.squaredStates      = 'even';
    pars_out.reservoirStateInit = 'random';
    pars_out.inputMatrixType    = 'balancedSparse';
    pars_out.inAmplitude        = 1.0;

    % overwrite pars_out values with pars_in values
    assert(isstruct(pars_in));
    names = fieldnames(pars_in);
    for k = 1:numel(names)
        pars_out.(names{k}) = pars_in.(names{k});
    end
end