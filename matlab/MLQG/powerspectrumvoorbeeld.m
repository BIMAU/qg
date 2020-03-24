nx = 128;
ny = 128;

xgrid = ((1:nx)-1)*2*pi/nx; % gridje zonder laatste punt
ygrid = ((1:ny)-1)*2*pi/ny;

E = sin(4*xgrid)'*sin(4*ygrid);
figure(1)
imagesc(E)

H = abs(fft2(E));
H = fftshift(H);

midx = nx/2+1;
midy = ny/2+1;

maxr = floor(sqrt(midx^2 + midy^2)) + 1;

rPrf = zeros(maxr,1);
cnt  = zeros(maxr,1);

C    = zeros(nx,ny);

for j = 1:ny
    for i = 1:nx
        r   = sqrt( (i-midx)^2 + (j-midy)^2 );
        idx = floor(r)+1;
        rPrf(idx) = rPrf(idx) + H(i,j);
        C(i,j) = idx;
        cnt(idx) = cnt(idx)+1;
    end
end

range = 2:maxr-4;
rPrf = rPrf(range);

figure(2)
loglog(range,rPrf,'k.-'); hold on;
%loglog(range,10*max(rPrf)*(range).^(-5/3));
%loglog(range,10*max(rPrf)*(range).^(-3));
hold off

vsm(C)
vsm(H)