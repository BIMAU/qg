nx=64
ny=nx
dest_time=0.9

stochsize=1000
  fprintf('IS YOUR STOCHSIZE CORRECT: %4d\n', stochsize)  
fid=fopen('timeslices.txt')
cnt=0;
ts=[];
%scroll through data file
while ((t <= dest_time ) &  ~feof(fid) )
  cnt=cnt+1;
  [V,Y,x,t,magnitudes]=rdtimeslice(fid,stochsize);
  t
end
  
fclose(fid);
[UY,SY,VY]=svd(Y',0);
%vsm(UY),vsm(VY),vsm(SY)
Var=(diag(SY)/sqrt(stochsize)).^2;
showV(V( :,1:4)*VY,Var(1:4),'DO',1,nx,nx,2,1,4)
pause

qg = QG(nx, ny);
%set Reynolds number
qg.set_par(5,10)
%set Full windstressfield
qg.set_par(11,1)
  
n = nx * ny * 2;

%Problem specification
%We want 
%We want a positive mass matrix since mass ought to be positive
M=qg.mass(n);
%Stochastic forcing
B=1e3*stochforcing(n,full(diag(Mass)));
A=qg.jacobian(x,0.0);


clear opts
opts.projection_method=1.2
opts.invA=@(x) A\x;
opts.restart_iterations=50;
%opts.restart_upon_convergence=false;
%opts.restart_tolerance=1e-6;
%opts.ortho='M'
%opts.restart=100
%opts.fast_orthogonalization=0
B=1e3*stochforcing(size(A,1),full(diag(M)));
size(B)
[S, MS, BS, Sinv, Vtrans] = RAILSschur(A, M, B);
size(BS)
fprintf('Reduced Lyapunov problem created\n')
opts.Ainv = Sinv;
[V,S,res,iter,resvec,timevec,restart_data] = RAILSsolver(S,MS,BS,1000,1e-13,opts);
plot(log10(resvec))
fprintf('Reduced Lyapunov system solved\n')
V = Vtrans(V);
if nx==32 
  X=V*S*V';
  norm(A*X*M'+M*X*A'+B*B')
end
fprintf('Basis prolongated\n')
%norm(V'*V-eye(size(V,2)))
%norm(V'*M*V-eye(size(V,2)))
  Var=diag(S(1:15,1:15))

showV(V( :,1:4),Var(1:4),'RAILS',1,nx,nx,2,1,4)

