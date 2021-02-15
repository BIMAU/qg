function [stopFlag, err, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)

    global memory

    use_fields   = false;
    use_spectrum = false;
    use_Ke       = true;
    % err_tol = 4.0;
    % err_tol = Inf;
    err_tol = 1.0e-4; % ~ 4*sigma

    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;

    windowsize = 10;

    if use_fields
        diff = abs(scaling*(testY-predY));
        err  = norm(diff(:),2) / norm(scaling*testY(:), 2);

    elseif use_spectrum
        testSpec = computeQGspectrum(qg, scaling*testY);
        predSpec = computeQGspectrum(qg, scaling*predY);

        sd   = find(testSpec > 1e-10);
        testSpec = testSpec(sd);
        predSpec = predSpec(sd);
        diff = (log(testSpec) - log(predSpec));
        err  = norm(diff,2);

    elseif use_Ke

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
        
        Ke_pred = sum(u2m_pred - um2_pred + v2m_pred - vm2_pred);
        Ke_test = sum(u2m_test - um2_test + v2m_test - vm2_test);
        err = abs(Ke_pred - Ke_test);
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