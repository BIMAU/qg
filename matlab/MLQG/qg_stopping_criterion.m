function [stopFlag, err, testSpec, predSpec] = qg_stopping_criterion(qg, predY, testY)
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;    
    
    testSpec = computeQGspectrum(qg, scaling*testY);
    predSpec = computeQGspectrum(qg, scaling*predY);

    sd   = find(testSpec > 1e-10);
    testSpec = testSpec(sd);
    predSpec = predSpec(sd);
    diff = (log(testSpec) - log(predSpec));
    err  = norm(diff,2);

    % ARBITRARY, make this more informed #FIXME
    if err > 4
        stopFlag = true;
    else
        stopFlag = false;
    end
end