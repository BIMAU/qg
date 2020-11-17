function [stopFlag, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;    
    
    testSpec = computeQGspectrum(qg, scaling*testY);
    predSpec = computeQGspectrum(qg, scaling*predY);

    sd   = find(testSpec > 1e-10);
    diff = (log(testSpec(sd)) - log(predSpec(sd)));
    err  = norm(diff,2);
    fprintf(' error: %1.2e\n', err);

    % ARBITRARY 
    if err > 5
        stopFlag = true;
    else
        stopFlag = false;
    end
end