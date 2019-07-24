function [Etdis, EtmeanV]=energy_transfers(V,x,vars,coskewness,Re)
[n,m]=size(V);
nx=sqrt(n/2)
dx=1/nx ;dy=dx;

%we need psi as a 2D field
psi=get_2Dpsifield(x,n,nx);
%psi contains zeros around the field
%compute staggered u,v
[u,v]=compute_uv(psi,nx,dx,dy);
%Compute grads at cell centers.
[muxc,muyc,mvxc,mvyc]=centergrads(u,v,nx,dx,dy);
%average velocities to cell centers
[u,v]=average_uv(u,v,nx);
uV(nx-1,nx-1,m)=0; vV=uV;
uxc(nx-1,nx-1,m)=0;uyc=uxc;vxc=uxc;vyc=uxc;
for k=1:m
  psi=get_2Dpsifield(V(:,k),n,nx);
  [uVtmp,vVtmp]=compute_uv(psi,nx,dx,dy);
  [uxc(:,:,k),uyc(:,:,k),vxc(:,:,k),vyc(:,:,k)]=centergrads(uVtmp,vVtmp,nx,dx,dy);
  [uV(:,:,k),vV(:,:,k)]=average_uv(uVtmp,vVtmp,nx);
  
end
for k=1:m
  Etdis(k)= -vars(k)/Re*dx*dy*sum(reshape( ...
    uxc(:,:,k).^2+uyc(:,:,k).^2+vxc(:,:,k).^2+vyc(:,:,k).^2,(nx-1)^2,1));   
  EtmeanV(k)= -vars(k)*dx*dy*sum(reshape( ...     
   muxc.*uV(:,:,k).^2+mvyc.*vV(:,:,k).^2+(muyc+mvxc).*vV(:,:,k).*uV(:,:,k),(nx-1)^2,1));      
end

end 
function [u,v]=compute_uv(psi,nx,dx,dy)
u(nx,nx+1)=0; v(nx+1,nx)=0;
rng1=1:nx-1;
u(:,rng1+1)=-(psi(:,rng1+1)-psi(:,rng1))/dy; %u(1,:)=u(nx,:)=0
v(rng1+1,:)=(psi(rng1+1,:)-psi(rng1,:))/dx;  %v(:,1)=v(:,nx)=0 
%extend u in y direction
%for a free slip boundary condition copy the current border -> u_y=0 at the boundary
%for a no-slip boundary copy minus the current border -> average u = 0
%at top and bottom we have slip
u(:,1)=u(:,2); u(:,nx+1)=u(:,nx);
%left and righ we have no slip
v(1,:)=-v(2,:); v(nx+1,:)=-v(nx,:);
end
function [uxc,uyc,vxc,vyc]=centergrads(u,v,nx,dx,dy)
%S1=u_x, S2=u_y+v_x, note that v_y=-u_
rng1=1:nx-1;
uxc=(u(rng1+1,rng1+1)-u(rng1,rng1+1))/dx; %no zeros around anymore
vyc=(v(rng1+1,rng1+1)-v(rng1+1,rng1))/dy; %no zeros around anymore
rng=1:nx;
tmp=(u(:,rng+1)-u(:,rng))/dy;
tmp=(tmp(:,rng1+1)+tmp(:,rng1))/2;
uyc=(tmp(rng1+1,:)+tmp(rng1,:))/2;
clear tmp
tmp=(v(rng+1,:)-v(rng,:))/dx;
tmp=(tmp(rng1+1,:)+tmp(rng1,:))/2;
vxc=(tmp(:,rng1+1)+tmp(:,rng1))/2;
end


function [u,v]=average_uv(u,v,nx)
rng1=1:nx-1;
u=(u(rng1+1,rng1+1)+u(rng1,rng1+1))/2; 
v=(v(rng1+1,rng1+1)+v(rng1+1,rng1))/2;  
end

function psi=get_2Dpsifield(x,n,nx)
for j=1:nx
  i=1:nx;
  psi(i,j)=x((j-1)*2*nx+2*i);
end
end 

