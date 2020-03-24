% Fine nx=256 data with fixed time step 
% data = load('N256_Re4.0e+04_Tstart141_Tend142_F0.5.mat');

T = size(data.states, 2);

for i = 1:T
    norm(data.states(:,i))
    data.times(i)
end



subplot(1,2,1)
plotQG(nx,ny,1,scaling*xf,false)
caxis(crange);
titleString = sprintf('Fine vorticity (day^{-1}), t = %3.0fd', ...
                      (data.times(idx)-data.times(1)) / day);
title(titleString);

subplot(1,2,2)
plotQG(nxc, nyc, 1, scaling*(F_a-F), false);
colorbar
caxis(Frange);
title('Eddy forcing for coarse model')
