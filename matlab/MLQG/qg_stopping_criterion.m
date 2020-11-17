function [stop] = qg_stopping_criterion(qg, predY, testY)
    Ldim    = 1e6;
    Udim    = 3.171e-2;
    tdim    = Ldim / Udim; % in seconds
    scaling = 3600*24/tdim;    
    
    testSpec = computeQGspectrum(qg, scaling*testY);
    fprintf('testSpec: %f\n', norm(testSpec));
end