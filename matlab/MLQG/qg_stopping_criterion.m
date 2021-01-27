function [stopFlag, err, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)
    use_fields   = false;
    use_spectrum = true;
    %err_tol = 4.0;
    err_tol = Inf;
    
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;    
    
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
    end

    if err > err_tol
        stopFlag = true;
    else
        stopFlag = false;
    end
end