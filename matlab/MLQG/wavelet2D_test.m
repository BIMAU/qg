% This is a demonstration of the nested permutations that can be used
% with a 2D block decomposition. Applying a wavelet to a column-major
% block permuted array gives non-local differences. This is visible in
% the stripe pattern (figure 3). With nested block permutations the
% differences are increasingly more local.

clear all
bs  = 32;       % block size for permutation P
n   = 128;      % grid size
nun = 2;        % number of unknowns in per grid point
dim = nun*n*n;  % total dimension

[all, QQ] = blockpermutation(n, n, nun, bs, true);
P = all{1};

% create data
A =  10*peaks(n/2);
A =  [A, A; A, A];
B = -peaks(n/2).^2/100;
B =  [B, B; B, B];
C =  peaks(n);

subplot(2,2,1)
imagesc(A')
title('dof 1')
set(gca,'ydir','normal')
colorbar

subplot(2,2,2)
imagesc(B')
title('dof 2')
set(gca,'ydir','normal')
colorbar

% get data in the ordering we're used to
if nun == 2
    A = [A(:), B(:)]';
elseif nun == 3
    A = [A(:), B(:), C(:)]';
end
A =  A(:);

% wavelet operator
H   = haarmat(bs^2);
H   = kron(speye(dim / (bs*bs)), H);

% reordering permutation
T  = speye(dim);
id = [];
for i = 1:bs^2
    id = [id, (i:bs^2:dim)];
end
T(:,id) = T(:,1:dim);

% dimension reduction
drd = dim/32;

% ordinary wavelet transform with single 'large' block permutation
Ht  = T*H*P;
% reduced dimension
Hp  = Ht(1:drd,:);
% compute reduced order coordinates / coefficients
x1  = Hp*A(:);
% transform back to original space
y1  = Hp'*x1;

subplot(2,2,3)
imagesc(reshape(y1(1:nun:dim),n,n)');
title('ordinary wavelet projection')
set(gca,'ydir','normal')
colorbar

% the same thing now with the full nested block permutation QQ
Ht  = T*H*QQ;
Hp  = Ht(1:drd,:);
x2  = Hp*A(:);
y2  = Hp'*x2;

subplot(2,2,4)
imagesc(reshape(y2(1:nun:dim),n,n)');
title('wavelet projection with nested P')
set(gca,'ydir','normal')
colorbar