function [V,Y,x,t,magnitudes]=rdtimeslice(fid,stochsize)
%
t   = fscanf(fid, '%e', 1);
%
n = fscanf(fid, '%d,%d',2);
m=n(2);
n=n(1);
%
magnitudes=fscanf(fid,'%e', m+1);
x = fscanf(fid, '%e', n);
V  = fscanf(fid, '%e', [n,m]);
Y = fscanf(fid, '%e', [m,stochsize]);
%
end
