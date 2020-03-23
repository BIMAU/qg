function [P] = blockpermutation(n,m,nun,bs)
    assert(mod(n,bs)==0)
    assert(mod(m,bs)==0)
    
    ndim = n*m*nun;
    P = sparse(ndim,ndim);
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
