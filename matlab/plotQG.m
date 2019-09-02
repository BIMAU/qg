function [] = plotQG(n,m,XX,state)
% Plotting QG model:
% n,m: horizontal grid size
% XX: 1,2 index of field you want to plot

% It's QG so nun = 2
nun = 2;
plotfield = reshape(state(XX:nun:end),n,m);

% constants 
udim = 1.6e-02; 
ldim = 1.0e+06; 
hdim = 6.0e+02; 
fact = udim*hdim*ldim/1.0e+06; 

% grid
for i=1:n
    x(i) = (i-1)/(n-1);
end

for j=1:m
    y(j) = (j-1)/(m-1); 
end

% scaling
maxp = max(max(plotfield)); 
% plot
contourf(x,y,fact*plotfield'/maxp); 
xlabel('x/L')
ylabel('y/L')
title('')

end