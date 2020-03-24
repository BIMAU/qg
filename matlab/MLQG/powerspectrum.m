
qg = QG(nx, ny, 1);

% QG params:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365 * day;

N = size(states,2)
for k = 1:ceil(N/10):N
    x = states(:,k);
    subplot(1,2,1)
    plotQG(nx,ny,1,3600*24/tdim*x,false);
    title(['relative vorticity \omega (day^{-1}), year: ', ...
           num2str(times(k)/year)])

    [u,v] = qg.compute_uv(x);

    u = reshape(u,nx,ny);
    v = reshape(v,nx,ny);
    om = reshape(x(1:2:end),nx,ny);
    ps = reshape(x(2:2:end),nx,ny);
       
    H = abs(fft2(u)).^2+abs(fft2(v)).^2;
    H = abs(fft2(u)).^2;
    
    %xgrid = ((1:nx)-1)*2*pi/nx;
    %ygrid = ((1:ny)-1)*2*pi/ny;
    %E = cos(5*xgrid)'*cos(5*ygrid);
    %plot(E(:,1),'k.-')
    %imagesc(E)
    %H = abs(fft2(E)).^2;
     
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
            r   = sqrt( (i-midx)^2 + (j-midy)^2 );
            idx = floor(r)+1;
            rPrf(idx) = rPrf(idx) + H(i,j);
            C(i,j) = idx;
            cnt(idx) = cnt(idx)+1;
        end
    end

    range = 2:maxr-3;
    %rPrf = rPrf./cnt;
    rPrf  = rPrf(range);
    
    % wavenumber related to different forcing types
    kx = 5;
    kx = 16;
    ky = kx;
    wavNum = floor(sqrt(kx^2+ky^2))+1;
    
    subplot(1,2,2)
    y1 = rPrf(wavNum-range(1)+1);
    y2 = rPrf(wavNum-range(1)+1);
    x1 = y1/(wavNum^-(5/3));
    x2 = y2/(wavNum^-(3));
    x3 = y2/(wavNum^-(5));
    loglog(range,rPrf,'k.-'); hold on;
    frst = 1:wavNum-range(1)+1;
    last = wavNum-range(1)+1:numel(range);
    loglog(range(frst),x1*(range(frst)).^(-5/3));
    loglog(range(last),x2*(range(last)).^(-3));
    %loglog(range(last),x3*(range(last)).^(-5));
    loglog([wavNum,wavNum],ylim,'k--');
    hold off
    xlim([min(range),max(range)]);
    input('');
end