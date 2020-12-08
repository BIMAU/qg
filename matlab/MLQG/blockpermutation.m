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

    if (nargin < 5) || ((bs / 2) < 2)
        nested = false;
    end

    assert(mod(n,bs)==0);
    assert(mod(m,bs)==0);

    ndom = n;
    mdom = m;
    sz = bs;
    P = create_permutation(ndom, mdom, nun, sz);
    Q = P;

    if (nested)
        % total number of block permutations
        total = 1;
        while ((sz / 2) >= 2)
            total = total + 1;
            sz  = sz / 2;
        end
        % if this fails the original bsz is not a power of 2
        assert(sz == 2);

        sz   = bs; % reset block sz
        out    = cell(total,1);
        out{1} = P;
        blocks = nun;
        for i = 2:total
            fprintf('Nested block permutation %d\n', i);

            % total number of blocks in full domain n*m
            blocks = blocks * (ndom / sz) * (mdom / sz);

            % domain in which we perform a new block decomposition
            ndom = sz; mdom = sz;

            % nested block sz
            sz = sz / 2;
            out{i} = kron(speye(blocks), ...
                          create_permutation(ndom, mdom, 1, sz));

            fprintf('  block sz: %d \n', sz);
            fprintf('  subdomain sz: %d x %d \n', ndom, mdom);
            fprintf('  number of blocks in whole domain: %d \n', blocks);

            Q = out{i} * Q; % combined operator
        end
        P = out;
    end
end

function [P] = create_permutation(n,m,nun,bs)
    dim = n*m*nun;
    P   = sparse(dim,dim);
    k   = 0;
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