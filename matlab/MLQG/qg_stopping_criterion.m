function [stopFlag, err, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)

    global memory windowsize
    
    use_fields   = false;
    use_spectrum = false;
    use_wavelet  = false;
    use_svd      = false;
    use_svdwav   = true;
    use_Ke       = false;
    use_Km       = false;
    use_Z        = false;
    use_mass     = false;
    use_ddtmass  = false;
    use_ddtPV    = false;

    err_tol = Inf;
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;

    nx = qg.nx;
    ny = qg.ny;

    region = [];
    for i = 0:6
        region = [region, i+(floor(nx/5)*nx:nx:1*nx/2*nx) - nx/2];
    end

    if use_fields
        % err_tol = 0.2;
        % diff = (testY-predY);
        % err  = norm(diff(:),2) / norm(testY(:), 2);
        
        err_tol  = 1;
                
        pred_coef = predY(:);
        test_coef = testY(:);

        memory = add_field(memory, 'pred_coef', pred_coef, windowsize);
        memory = add_field(memory, 'test_coef', test_coef, windowsize);

        T = size(memory.pred_coef,2);       
        
        if T > 1
            MSE   = mean(sum((memory.pred_coef - memory.test_coef).^2));
            NRM   = mean(sum((memory.test_coef - mean(memory.test_coef,2)).^2));
            NRMSE = sqrt(MSE / NRM);
            err   = NRMSE;
        else
            err = 0;
        end


    elseif use_spectrum
        err_tol  = 1.5;
        testSpec = computeQGspectrum(qg, testY);
        predSpec = computeQGspectrum(qg, predY);

        sd = 2:40;
        % pred_coef = log(predSpec(sd));
        % test_coef = log(testSpec(sd));

        pred_coef = (predSpec(sd));
        test_coef = (testSpec(sd));

        memory = add_field(memory, 'pred_coef', pred_coef, windowsize);
        memory = add_field(memory, 'test_coef', test_coef, windowsize);

        T = size(memory.pred_coef,2);

        if T > 1

            MSE   = mean(sum((memory.pred_coef - memory.test_coef).^2));
            NRM   = mean(sum((memory.test_coef - mean(memory.test_coef,2)).^2));
            NRMSE = sqrt(MSE / NRM);
            err   = NRMSE;
        else
            err = 0;
        end

    elseif use_wavelet
        err_tol  = 0.5;

        % pred_coef = log(abs(memory.Ha' * predY(:)));
        % test_coef = log(abs(memory.Ha' * testY(:)));
        pred_coef = ((memory.Ha' * predY(:)));
        test_coef = ((memory.Ha' * testY(:)));

        memory = add_field(memory, 'pred_coef', pred_coef, windowsize);
        memory = add_field(memory, 'test_coef', test_coef, windowsize);

        % diff = abs(sum(abs(pred_coef))-sum(abs(test_coef))) / abs(test_coef(1));
        % diff = sum(abs(pred_coef - test_coef)) / sum(abs(test_coef(:)));
        % diff = norm(pred_coef - test_coef) ./ abs(test_coef(1));
        % diff = abs(test_coef).*(pred_coef-test_coef).^2 / abs(test_coef(1));
        % diff = abs(sum( (test_coef).^2 ) - sum( (pred_coef).^2)) ./ sum( (test_coef).^2);
        
        T = size(memory.pred_coef,2);

        if T > 1
            MSE   = mean(sum((memory.pred_coef - memory.test_coef).^2));
            NRM   = mean(sum((memory.test_coef - mean(memory.test_coef,2)).^2));
            NRMSE = sqrt(MSE / NRM);
            err   = NRMSE;
        else
            err = 0;
        end

        % if T > 0
        %     ERR   = (pred_coef - test_coef).^2;
        %     VAR   = memory.VAR_wav;
        %     err   = sqrt(mean(ERR ./ VAR));
        % else
        %     err = 0;
        % end

    elseif use_svd
        err_tol  = 0.5;

        pred_coef = memory.U' * predY(:);
        test_coef = memory.U' * testY(:);

        memory = add_field(memory, 'pred_coef', pred_coef, windowsize);
        memory = add_field(memory, 'test_coef', test_coef, windowsize);

        % diff = abs(sum(abs(pred_coef))-sum(abs(test_coef))) / abs(test_coef(1));
        % diff = sum(abs(pred_coef - test_coef)) / sum(abs(test_coef(:)));
        % diff = norm(pred_coef - test_coef) ./ abs(test_coef(1));
        % diff = abs(test_coef).*(pred_coef-test_coef).^2 / abs(test_coef(1));
        % diff = abs(sum( (test_coef).^2 ) - sum( (pred_coef).^2)) ./ sum( (test_coef).^2);

        T = size(memory.pred_coef,2);
        
        % if T > 1

        %     pred_varcoef = pred_coef - mean(memory.pred_coef,2);
        %     test_varcoef = test_coef - mean(memory.test_coef,2);            
        %     memory = add_field(memory, 'pred_varcoef', pred_varcoef, windowsize);
        %     memory = add_field(memory, 'test_varcoef', test_varcoef, windowsize);
        %     MSE =  sum( ( ( mean(memory.pred_varcoef, 2).^2 ) - ...
        %                   ( mean(memory.test_varcoef, 2).^2 ) ).^2 ) ;
            
        %     err = MSE;
        % else
            
        %     err = 0;
            
        % end

        if T > 1
            MSE   = mean(sum((memory.pred_coef - memory.test_coef).^2));
            NRM   = mean(sum((memory.test_coef - mean(memory.test_coef,2)).^2));
            NRMSE = sqrt(MSE / NRM);
            err   = NRMSE;
        else
            err = 0;
        end

        % if T > 0
        %     ERR   = (pred_coef - test_coef).^2;
        %     VAR   = memory.VAR_svd;
        %     err   = sqrt(mean(ERR ./ VAR));
        % else
        %     err = 0;
        % end
    elseif use_svdwav
        err_tol  = 0.5;

        pred_coef = memory.Uwav' * memory.Ha' * predY(:);
        test_coef = memory.Uwav' * memory.Ha' * testY(:);

        memory = add_field(memory, 'pred_coef', pred_coef, windowsize);
        memory = add_field(memory, 'test_coef', test_coef, windowsize);

        T = size(memory.pred_coef,2);
        
        % figure(3);
        % for i = 1:16                                                         
        %     f = memory.Ha*memory.Uwav(:,i);                                      
        %     subplot(4,4,i);
        %     imagesc(reshape(f(1:2:end),64,64)'); set(gca,'ydir','normal')    
        % end
        
        % keyboard
        
        if T > 1
            MSE   = mean(sum((memory.pred_coef - memory.test_coef).^2));
            NRM   = mean(sum((memory.test_coef - mean(memory.test_coef,2)).^2));
            NRMSE = sqrt(MSE / NRM);
            err   = NRMSE;
        else
            err = 0;
        end

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
            err = abs(Ke_pred - Ke_test) / abs(Ke_test);

        elseif use_Km
            err_tol = 0.2;
            Km_pred = sum(um2_pred(region) + vm2_pred(region));
            Km_test = sum(um2_test(region) + vm2_test(region));
            err = abs(Km_pred - Km_test) / abs(Km_test);
            % err = 0; % never stop
        end

    elseif use_Z
        err_tol = 0.2;
        om_pred = predY(1:2:end)';
        om_true = testY(1:2:end)';
        Z_pred = sum(om_pred.^2);
        Z_true = sum(om_true.^2);
        err = abs(Z_pred - Z_true) / abs(Z_true);

    elseif use_mass
        err_tol = Inf; % never stop

        ps_pred = predY(2:2:end)';
        ps_true = testY(2:2:end)';
        M_pred  = sum(ps_pred);
        M_true  = sum(ps_true);

        err = abs(M_pred - M_true) / abs(M_true);

    elseif use_ddtmass
        err_tol = 2;

        ps_pred = predY(2:2:end)';
        ps_true = testY(2:2:end)';
        M_pred  = sum(ps_pred);
        M_true  = sum(ps_true);

        memory = add_field(memory, 'M_pred', M_pred, 2);
        memory = add_field(memory, 'M_true', M_true, 2);

        err = 0;
        if numel(memory.M_pred) == 2
            ddtM_pred = memory.M_pred(2) - memory.M_pred(1);
            ddtM_true = memory.M_true(2) - memory.M_true(1);
            err = abs(ddtM_true - ddtM_pred);
        end

    elseif use_ddtPV

        err_tol = 10;

        beta = 504.5727;
        dt = 0.0027;

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

            err = abs(ddtPV_true-ddtPV_pred) / abs(ddtPV_true);
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