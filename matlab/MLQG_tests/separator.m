function [P, nDetails, nAverage] = separator(ndim, blocksize, cutoff)

% create separation permutation

    P = speye(ndim);
    ord1   = [];
    ord2   = [];
    b2 = blocksize^2;
    for i = 1:ndim/b2
        ord1 = [ord1, (i-1)*b2+1:(i-1)*b2+cutoff];
        ord2 = [ord2, (i-1)*b2+cutoff+1:i*b2];
    end

    P = P([ord1,ord2],:);
    nAverage = ndim/b2*cutoff;
    nDetails = ndim - nAverage;
end