addpath /mnt/D/RAILS/matlab
addpath /home/p167800/swe_perturbed/matlab
addpath /project/fwubs/fredwubs/matlab
var = mmread('VarFile.mm');
V = mmread('V.mm');
nx=64;
DS=diag(var);
if 1
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
end
if 0
%Morthogonalize V
[V,R]=qrM(V,M);
[P,Vd]=eig(R*DS*R');
[Var,iVar]=sort(diag(Vd),'descend');
fprintf('Variances for M-orthogonal basis')
Var'
V=V*P(:,iVar);
end
m=4
showV(V(:,1:m),Var(1:m),'DO',1,nx,nx,2,1,m)
print -bestfit -dpdf streamFlds.pdf
figure(4)
print -bestfit -dpdf vorticityFlds.pdf
quit
