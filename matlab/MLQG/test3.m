clear all

IM    = flipud(load_image())';
nun   = 2;
[n,m] = size(IM);
ndim  = nun*m*n;

mxIm  = max(IM(:));

x = zeros(ndim,1);
ord = [];

figure(1);
imagesc(IM');
set(gca,'ydir','normal');


for j = 1:n
    for i = 1:m
        ord = [ord, nun*(m*(j-1)+(i-1))+1];
        
        x(ord(end)) = IM(i,j);
        
        ord = [ord, nun*(m*(j-1)+(i-1))+2];
        
        x(ord(end)) = mxIm-IM(i,j);
    end
end
        
bs = 4; % blocksize

IM = reshape(x, nun, m, n);
colormap(gray);
figure(1);
imagesc(squeeze(IM(1,:,:))'); set(gca,'ydir','normal');
figure(2);
imagesc(squeeze(IM(2,:,:))'); set(gca,'ydir','normal');
colormap(gray);

%---------------------------------------------

P1 = blockpermutation(n,m,nun,bs);
z  = P1*x;

IM = reshape(z(1:bs^2),bs,bs);
figure(3)
imagesc(IM'); set(gca,'ydir','normal');
colormap(gray);

IM = reshape(z(bs^2+1:2*bs^2),bs,bs);
figure(4)
imagesc(IM'); set(gca,'ydir','normal');
colormap(gray);

H = haarmat(bs^2);
M = speye(ndim/(bs^2));
H = kron(M,H);
w = H*z;

[P2, nDetails, nAverage] = separator(ndim, bs, 1);

v = P2*w;
average = v(1:nAverage);
details = v(nAverage+1:end);
details(abs(details)<40) = 0;
v = [average; details];
Q = P2*H*P1;
X = Q'*v;

fprintf('compression: %f\n', ndim/nnz(v));

IM = reshape(X, nun, m, n);
colormap(gray);
figure(5);
imagesc(squeeze(IM(1,:,:))'); set(gca,'ydir','normal');
colormap(gray);
figure(6);
imagesc(squeeze(IM(2,:,:))'); set(gca,'ydir','normal');
colormap(gray);

function [IM] = load_image()
% load image
    I       = imread('frank.jpg');
    I       = mean(I,3);
    factor  = 4;
    I       = I(1:factor:end,1:factor:end);
    [m,n]   = size(I);
    x       = linspace(0,1,n); y = linspace(0,1,m);
    [X,Y]   = meshgrid(x,y);
    t       = nextpow2(max([m,n]));
    xi      = linspace(0,1,2^t); yi = linspace(0,1,2^t);
    [Xi,Yi] = meshgrid(xi,yi);
    IM      = interp2(X,Y,I,Xi,Yi);
end

function [bool] = in_range(i,j,x,rangei,rangej,rangex)    
    bool = false;
    bool = ((i >= min(rangei)) && ...
            (i <= max(rangei)) && ...
            (j >= min(rangej)) && ...
            (j <= max(rangej)) && ...
            (x >= min(rangex)) && ...
            (x <= max(rangex)))    
end