;
nx=64

stochsize=1000
%stochsize=2
  fprintf('IS YOUR STOCHSIZE CORRECT: %4d\n', stochsize)  
%fid=fopen('timeslicesOddEvenDwsymStochForc01.txt')
fid=fopen('timeslices.txt')

while ~feof(fid)
  [V,Y,x,t,magnitudes]=rdtimeslice(fid,stochsize);
  t
  size(x)
    [n,m]=size(V)
  size(Y)
  %Y=Y';
  %pause
  excess=0; 
  %Vt=comp_moments(Y,4,excess);
  excess=1
  [Vt,expect,moments,Var]=comp_moments(Y,4,excess);
  Y=Vt'*Y;
  CoSkewn=1; CoKurt=1; CoExKurt=moments.m4;
  %plotdistr(Y,expect,Var,CoSkewn, CoKurt, CoExKurt)

  V=V*Vt;

    DS=diag(Var);
   
    D=spdiags(kron([-1 2 -1],ones(nx,1)),[-1,0,1],nx,nx);
   %Since psi has zeroes all around the first and last column/row should be made zero 
   D([1,nx],:)=0;D(:,[1,nx])=0;
   D=kron(speye(nx),D)+kron(D,speye(nx)); D=kron(D,sparse([0 0;0 1]));
   [V,R]=qrM(V,D);
   [P,Vd]=eig(R*DS*R');
   [Var,iVar]=sort(diag(Vd),'descend');
   fprintf('Variances for Laplacian-orthogonal basis\n')
  Var'
   %(R*DS*R')P=P*Vd -> V*R*DS*R'*V'= V*P*Vd*P'V' 
   V=V*P(:,iVar);
  m=2
  showV([x,V],[1;Var],'DO-mode',1,nx,nx,2,1,m)
  %pause
  %showV(x,1,'Average sol',1,nx,nx,2,1,1)
  %pause
  Y=P(:,iVar)'*(R*Y);
    
  plotdistr(Y,expect,Var,CoSkewn, CoKurt, CoExKurt)
  pause
  
end
fclose(fid);


%energies(V,sol,vars,coskewness)
