function [stopFlag, err, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)

    global memory windowsize

    use_fields   = true;
    use_spectrum = false;
    use_wavelet  = false;
    use_svd      = false;
    use_svdwav   = false;
    use_Ke       = false;
    use_Km       = false;
    use_Z        = false;
    use_mass     = false;
    use_ddtmass  = false;
    use_ddtPV    = false;

    err     = 0;
    err_tol = Inf;
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;

    beta = 504.5727;
    dt   = 0.0027;

    nx = qg.nx;
    ny = qg.ny;

    region = [];
    for i = 0:6
        region = [region, i+(floor(nx/5)*nx:nx:1*nx/2*nx) - nx/2];
    end

    if use_fields
        err_tol  = 1.0;
        err = NRMSE(predY(:), testY(:));
        
    elseif use_spectrum
        err_tol  = 1.0;
        testSpec = computeQGspectrum(qg, testY);
        predSpec = computeQGspectrum(qg, predY);
        
        % numel(testSpec)
        sd = 2:numel(testSpec)-5;
         pred_coef = log(predSpec(sd));
         test_coef = log(testSpec(sd));
         % pred_coef = (predSpec(sd));
         % test_coef = (testSpec(sd));

        err = NRMSE(pred_coef, test_coef);

    elseif use_wavelet
        err_tol  = 1.0;

        %pred_coef = log(abs(memory.Ha' * predY(:)));
        %test_coef = log(abs(memory.Ha' * testY(:)));
        pred_coef = memory.Ha' * predY(:);
        test_coef = memory.Ha' * testY(:);

        err = NRMSE(pred_coef, test_coef);

    elseif use_svd
        err_tol  = 1.0;

        pred_coef = memory.U' * predY(:);
        test_coef = memory.U' * testY(:);

        err = NRMSE(pred_coef, test_coef);

    elseif use_svdwav
        err_tol  = 1.0;

        pred_coef = memory.Uwav' * memory.Ha' * predY(:);
        test_coef = memory.Uwav' * memory.Ha' * testY(:);

        err = NRMSE(pred_coef, test_coef);

    elseif use_Ke || use_Km

        [u_pred, v_pred] = qg.compute_uv(scaling*predY);
        [u_test, v_test] = qg.compute_uv(scaling*testY);

        memory = add_field(memory, 'u_pred', u_pred, windowsize);
        memory = add_field(memory, 'v_pred', v_pred, windowsize);
        memory = add_field(memory, 'u_test', u_test, windowsize);
        memory = add_field(memory, 'v_test', v_test, windowsize);

        u2m_pred = mean(memory.u_pred.^2, 2);
        um2_pred = mean(memory.u_pred, 2).^2;
        v2m_pred = mean(memory.v_pred.^2, 2);
        vm2_pred = mean(memory.v_pred, 2).^2;
        u2m_test = mean(memory.u_test.^2, 2);
        um2_test = mean(memory.u_test, 2).^2;
        v2m_test = mean(memory.v_test.^2, 2);
        vm2_test = mean(memory.v_test, 2).^2;

        if use_Ke
            err_tol = 0.5;

            Ke_pred = sum(u2m_pred(region) - um2_pred(region) + v2m_pred(region) - vm2_pred(region));
            Ke_test = sum(u2m_test(region) - um2_test(region) + v2m_test(region) - vm2_test(region));

            err = NRMSE(Ke_pred, Ke_test);

        elseif use_Km
            err_tol = 0.5;

            Km_pred = sum(um2_pred(region) + vm2_pred(region));
            Km_test = sum(um2_test(region) + vm2_test(region));

            err = NRMSE(Km_pred, Km_test);
        end

    elseif use_Z

        err_tol = 2.0;
        om_pred = predY(1:2:end)';
        om_true = testY(1:2:end)';
        Z_pred = sum(om_pred.^2);
        Z_true = sum(om_true.^2);

        err = NRMSE(Z_pred, Z_true);

    elseif use_mass
        err_tol = 0.5;

        ps_pred = predY(2:2:end)';
        ps_true = testY(2:2:end)';
        M_pred  = sum(ps_pred);
        M_true  = sum(ps_true);

        err = NRMSE(M_pred, M_true);

    elseif use_ddtmass
        err_tol = 0.5;

        ps_pred = predY(2:2:end)';
        ps_true = testY(2:2:end)';
        M_pred  = sum(ps_pred);
        M_true  = sum(ps_true);

        memory = add_field(memory, 'M_pred', M_pred, 2);
        memory = add_field(memory, 'M_true', M_true, 2);

        err = 0;
        if numel(memory.M_pred) == 2
            ddtM_pred = (memory.M_pred(2) - memory.M_pred(1))/dt;
            ddtM_true = (memory.M_true(2) - memory.M_true(1))/dt;

            err = NRMSE(ddtM_true, ddtM_pred);
        end

    elseif use_ddtPV
        
        err_tol = 0.5;

        om_true = testY(1:2:end);
        om_pred = predY(1:2:end);

        memory = add_field(memory, 'om_true', om_true, 2);
        memory = add_field(memory, 'om_pred', om_pred, 2);

        err = 0;

        if size(memory.om_true,2) == 2

            ddtom_true = (memory.om_true(:,2) - memory.om_true(:,1)) / dt;
            ddtom_pred = (memory.om_pred(:,2) - memory.om_pred(:,1)) / dt;

            [omx_true, omy_true] = qg.gradient(om_true(:));
            [u_true, v_true]     = qg.compute_uv(testY);
            [omx_pred, omy_pred] = qg.gradient(om_pred(:));
            [u_pred, v_pred]     = qg.compute_uv(predY);

            advom_true = u_true.*omx_true + ...
                v_true.*omy_true;
            advom_pred = u_pred.*omx_pred + ...
                v_pred.*omy_pred;

            bv_true = beta*v_true;
            bv_pred = beta*v_pred;
            
            ddtPV_true_full = ddtom_true + advom_true + bv_true;
            ddtPV_pred_full = ddtom_pred + advom_pred + bv_pred;
            
            ddtPV_true = sum(ddtPV_true_full(region));
            ddtPV_pred = sum(ddtPV_pred_full(region));

            err = NRMSE(ddtPV_true, ddtPV_pred);
        end

    end
    
    if err > err_tol
        stopFlag = true;
    else
        stopFlag = false;
    end
end


function [mem] = add_field(mem, name, field, wsize)
    if isfield(mem, name)
        % append to mem
        mem.(name) = [mem.(name), field];

        % retain max wsize fields
        memsize = size(mem.(name), 2);
        mem.(name) = mem.(name)(:, max(1,memsize-wsize+1):memsize);
    else
        % create field
        mem.(name) = field;
    end
end

function [err, NRM] = NRMSE(pred, test)

    global memory windowsize

    memory = add_field(memory, 'pred', pred, windowsize);
    memory = add_field(memory, 'test', test, windowsize);

    diff =  memory.pred - memory.test;
    
    tvar =  memory.test - mean(memory.test, 2);

    T   = size(diff,2);
    N   = size(diff,1);
    NRM = 1;
    
    if T < windowsize        
        % padding diff with zeros
        diff = [zeros(N, windowsize-T), diff];
    end
    
    if T > 1
        MSE   = mean(sum( diff.^2 , 1));
        NRM   = mean(sum( tvar.^2 , 1));
        NRMSE = sqrt(MSE / NRM);
        err   = NRMSE;
    else
        err = 0;
    end

end