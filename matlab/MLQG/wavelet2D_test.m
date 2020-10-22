% This is a demonstration of the nested permutations that can be used
% with a 2D block decomposition. Applying a wavelet to a column-major
% block permuted array gives non-local differences. This is visible in
% the stripe pattern (figure 3). With nested block permutations the
% differences are increasingly more local.

bs = 32;        % block size for permutation P
n   = 128;      % grid size
nun = 2;        % number of unknowns in per grid point
dim = nun*n*n;  % total dimension

[all, QQ] = blockpermutation(n, n, nun, bs, true);

% create data
A =  peaks(n);
B = -peaks(n).^2/100;

figure(1); 
imagesc(A')
set(gca,'ydir','normal')
colorbar

figure(2); 
imagesc(B')
set(gca,'ydir','normal')
colorbar

% get data in the ordering we're used to
A = [A(:), B(:)]';
A =  A(:);

% wavelet operator
H = haarmat(bs^2);
H = kron(speye(dim / bs / bs), H);

nested = size(all,1);

r = randn(dim,1);
PP = speye(dim);
for i = 1:nested
    PP   = all{i} * PP ;
    Hrqp = T*H*PP;
    figure(2+i)
    x = Hrqp*A(:);
    x = x + r.*x / 10;
    y = Hrqp'*x;
    plotQG(n,n,2,y,false)
end