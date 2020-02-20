function D=innerprodmatrix(nx)
D=spdiags(kron([-1 2 -1],ones(nx,1)),[-1,0,1],nx,nx);
%Since psi has zeroes all around the first and last column/row should be made zero 
D([1,nx],:)=0;D(:,[1,nx])=0;
D=kron(speye(nx),D)+kron(D,speye(nx)); 
% Only apply it to the psi (order is zeta,psi,zeta,psi, ....)
D=kron(D,sparse([0 0;0 1]));
D=[]; 
end
