function [P,Q] = blockpermutation(n,m,nun,bs,nested)
% Block permute a column-major 2D array with several unknowns per
% node.
%
% input:
%   n:      number of columns
%   m:      number of rows
%   nun:    number of unknowns
%   bs:     block size (square)
%   nested: compute nested block permutations
%
% output:
%   P:      (cell array of) permutation operator(s)
%   Q:      combined permutation operator

    if nargin < 5
        nested = false;
    end

    assert(mod(n,bs)==0);
    assert(mod(m,bs)==0);

    ndom = n;
    mdom = m;
    size = bs;
    P = create_permutation(ndom, mdom, nun, size);
    Q = P;

    if (nested)
        % total number of block permutations
        total = 1;
        while ((size / 2) >= 2)
            total = total + 1;
            size = size / 2;
        end
        % if this fails the original bsize is not a power of 2
        assert(size == 2);

        size   = bs; % reset block size
        out    = cell(total,1);
        out{1} = P;
        blocks = 1;
        for i = 2:total
            fprintf('\nnested block permutation\n');

            % total number of blocks in full domain n*m
            blocks = blocks * (ndom / size) * (mdom / size);

            % domain in which we perform a new block decomposition
            ndom = size; mdom = size;

            % nested block size
            size = size / 2;
            out{i} = kron(speye(blocks), ...
                          create_permutation(ndom, mdom, nun, size));

            fprintf('block  size: %d \n', size);
            fprintf('domain size: %d x %d \n', ndom, mdom);
            fprintf('total number of blocks: %d \n', blocks);

            Q = out{i} * Q; % combined operator
        end
        P = out;
    end
end

function [P] = create_permutation(n,m,nun,bs)
    dim = n*m*nun;
    P = sparse(dim,dim);
    k = 0;
    for posj = 0:bs:n-bs
        rangej = posj+1:posj+bs;
        for posi = 0:bs:m-bs
            rangei = posi+1:posi+bs;
            for xx = 1:nun
                for j = rangej
                    for i = rangei
                        k = k + 1;
                        col = nun*(n*(j-1)+(i-1))+xx;
                        P(k,col) = 1;
                    end
                end
            end
        end
    end
    P = sparse(P);
end