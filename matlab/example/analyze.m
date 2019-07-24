nx=64
Re=40
stochsize=1000
%stochsize=2
  fprintf('IS YOUR STOCHSIZE CORRECT: %4d\n', stochsize)  
%fid=fopen('timeslicesOddEvenDwsymStochForc01.txt')
fid=fopen('timeslices.txt')
%  fid=fopen('timeslicesDt1e-3ndtsub10m4Frc1.txt')
%  fid=fopen('timeslicesDt1e-3ndtsub10m4.txt')
%  fid=fopen('timeslicesSymForcing_dt1e-3_ndtsub10_m4.txt')
%timeslicesDt1e-3ndtsub10m4Frc10.txt')
  cnt=0;
ts=[];
while ~feof(fid)
  close all
  cnt=cnt+1;
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
  m=4
  Y=P(:,iVar)'*(R*Y);
  showV([x,V],[0;Var],'Solution and DO-fields',1,nx,nx,2,1,m)
  %pause
  indYp=find(Y(1,:) > 1.0); length(indYp)
  indYn=find(Y(1,:) < -1.0); length(indYn)
  Yp=Y(:,indYp); expectp=sum(Yp,2)/length(indYp)  
  Yn=Y(:,indYn); expectn=sum(Yn,2)/length(indYn) 
 % showV([x+expectp(1)*V(:,1),x+expectn(1)*V(:,1)],[0;0],'Approx. Stable solutions',1,nx,nx,2,1,2)
  %showV(x,1,'Average sol',1,nx,nx,2,1,1)
  %pause
  %Y=P(:,iVar)'*(R*Y);
  
  excess=0; 
  %Vt=comp_moments(Y,4,excess);
  excess=1
  [Vt,expect,moments,Var]=comp_moments(Y,4,excess);
  Y=Vt'*Y;
  V=V*Vt;
  Vt % should be identity
  
  Var
  CoSkewn=moments.m3 
  CoKurt=1; 
  CoExKurt=moments.m4
  mV=size(V,2);
  [Etdis(1:mV,cnt),EtmeanV(1:mV,cnt)]=energies(V,x,Var,CoSkewn,Re) 
  ts=[ts,t];
  plotdistr(Y,expect,Var,CoSkewn, CoKurt, CoExKurt)
%  if (t==7) 
    pause
%  end
  
end
fclose(fid);



