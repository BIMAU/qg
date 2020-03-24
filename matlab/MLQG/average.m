function [out] = average(in, nx, ny, nun, ff)
    assert(mod(nx*ny*nun, ff) == 0);
    
    in = reshape(in,nun,nx,ny);


    
    nxc = nx / ff;
    nyc = ny / ff;

    out = zeros(nun, nxc, nyc);
    
    idj = 1:nyc;
    idi = 1:nxc;
    
    %out(1,idi,idj) = sum(sum(in(1, ff*(idi-1)+(1:ff), ff*(idj-1)+(1:ff))))/(ff*ff);
    %out(2,idi,idj) = sum(sum(in(2, ff*(idi-1)+(1:ff), ff*(idj-1)+(1:ff))))/(ff*ff);
    keyboard
    
    for j = 1:nyc
        for i = 1:nxc
            out(1,i,j) = sum(sum(in(1, ff*(i-1)+(1:ff), ff*(j-1)+(1:ff))))/(ff*ff);
            out(2,i,j) = sum(sum(in(2, ff*(i-1)+(1:ff), ff*(j-1)+(1:ff))))/(ff*ff);
        end
    end
    
    out = out(:);
end