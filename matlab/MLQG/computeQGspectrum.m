function [rPrf, C, maxr, sumCoefs] = computeQGspectrum(qg, x)
    
    nx = qg.nx;
    ny = qg.ny;
    
    assert(numel(x) == nx*ny*qg.nun);
    
    % get velocity field    
    [u,v] = qg.compute_uv(x(:));
    Udim  = qg.Udim;
    u     = Udim*reshape(u,nx,ny);
    v     = Udim*reshape(v,nx,ny);
    
    H = abs(fft2(u)).^2+abs(fft2(v)).^2;
    
    H = fftshift(H);
    rows = size(H,1);
    cols = size(H,2);
    
    midx = rows/2+1;
    midy = cols/2+1;

    maxr = floor(sqrt(midx^2 + midy^2)) + 1;

    rPrf = zeros(maxr,1);
    cnt  = zeros(maxr,1);
    C    = zeros(rows,cols);

    for j = 1:cols
        for i = 1:rows
            r         = sqrt( (i-midx)^2 + (j-midy)^2 );
            idx       = floor(r)+1;
            rPrf(idx) = rPrf(idx) + H(i,j);
            C(i,j)    = idx;
            cnt(idx)  = cnt(idx)+1;
        end
    end
    
    sumCoefs = sum(rPrf);
end