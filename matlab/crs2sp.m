function [A] = crs2sp(beg,jco,co)

% Create sparse matrix based on 1-based CRS arrays

    n   = numel(beg)-1;
    nnz = beg(end)-1;

    ivals = zeros(nnz,1);
    jvals = jco;
    vals  = co;

    row = 1;
    idx = 1;
    while row <= n
        for k = beg(row):beg(row+1)-1
            ivals(idx) = row;
            idx        = idx + 1;
        end
        row = row + 1;
    end
    A = sparse(ivals, jvals, vals, n, n);
end

