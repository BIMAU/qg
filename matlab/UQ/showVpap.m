function []=showV(V,var,method,nz,nx,ny,blksize,CASE,mp,fig_offset);
% []=showV(V,var,method,nz,nx,ny,blksize,CASE,mp);
% V  : nx*nz*ny*blksize x m matrix;
% var: m variances for printing above the plot
% method : string to be added on top of the plot
% nz : number of layers
% nx : number of points in x-direction
% ny : number of points in y-direction
% blksize : number of unknowns per cell/gridpoint
% CASE : indicates which crosscut plane should be visualized
% mp: number of vectors in V to be visualized.

n=size(V,1);
m=size(V,2);
n=nx*nz*ny*blksize; 
if ( m!= 4); return 'm neq 4'
[ha,hfig]=tight_subplot_cm(2,2,[2.0,1.9],[2.8,1.4],[1.8,.8],18,25,'No');
for i=1:blksize
  f=figure(fig_offset+i+3);
  f.Units='centimeters';
  f.OuterPosition=[0,0,22,29];
  f.PaperPositionMode='auto';
  f.PaperType='A5';
  p = uipanel('Parent',f,'BorderType','none'); 
  p.Title = [method, sprintf(' ; nx=ny=  %d',nx)]; 
  p.TitlePosition = 'centertop'; 
  p.FontSize = 12;
  p.FontWeight = 'bold';
  for j=1:mp
    fld=reshape(V(i:blksize:n,j),nx,ny,nz);
    switch CASE 
    case 1
      restruct=2;
      if (mp==1) restruct=1; end;
      for k=1:nz    
        subplot(ceil(mp/restruct),restruct,k+(j-1)*nz,'Parent',p);
        contourf(transpose(fliplr(fld(:,:,k))),40)
        axis equal
	title(sprintf('Variance %1.5g', var(j)))
	colorbar
      end
    case 2
      nxl=1
      for k=1:nxl
	subplot(mp,nxl,k+(j-1)*nxl);
	fldtmp(:,:)=fld(k,:,:);
	contourf(fldtmp')
      end
    case 3
      for k=1:ny
	subplot(mp,ny,k+(j-1)*ny);
	fldtmp(:,:)=fld(:,k,:);
	contourf(fldtmp)
      end
    end
  end
end
