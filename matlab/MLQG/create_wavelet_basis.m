function [H_out] = create_wavelet_basis(nx, ny, nun, bs, nested)
% Create wavelet basis for a 2D grid with nested block decompositions
% #TODO make nesting level optional

    time = tic;
    fprintf('Create wavelet basis... \n');

    dim = nx*ny*nun;

    [all, P] = blockpermutation(nx,ny,nun,bs,nested);

    H = haarmat(bs^2);

    % create a block diagonal matrix with dim/(bs^2) H blocks
    M = speye(dim / (bs^2));
    H = kron(M, H);

    % create row ordering: orders the average and difference modes from
    % coarse to fine
    T  = speye(dim);
    id = [];
    for i = 1:bs^2
        id = [id, (i:bs^2:dim)];
    end
    T(:,id) = T(:,1:dim);

    H_out = T*H*P;
    fprintf('Create wavelet basis... done (%f)\n', toc(time));
end