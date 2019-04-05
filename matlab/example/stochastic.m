nx = 40;
ny = 40;
n = nx * ny * 2;
qg = QG(nx, ny);

% Set zeta
zeta = 1/n;
qg.set_par(20, zeta);

%set Reynolds number
qg.set_par(5,10)
%set Full windstressfield
qg.set_par(11,1)
  
x = zeros(n, 1);
tol=1e-10;

%x=newtonfw(@(y)-qg.rhs(y),@(y)qg.jacobian(y,0.0),x,tol);
m=2;
stochiter=1000;
V=zeros(n,m);
Y=zeros(m,stochiter);
Mass=qg.mass(n);

integr_pars.ndtsub=1;
integr_pars.dt=1e-3;
integr_pars.T=1e-1;

c=1e5*stochforcing(n,full(diag(Mass)));

%c=1e-5*ones(n,1);
f1=@(X,b)b-qg.rhs(X);
f2=@(X)qg.jacobian(X,0.0);
f3=@(X,Y)qg.bilin(X,Y); 
[x,V,Y,timings,timeseries]=DONonlin(x,V,Y,f1,f2,f3,c,Mass,integr_pars);
timings
figure(1)
hold off
semilogy(timeseries.time,timeseries.nrmdetsol);
hold on
for i=1:m
  semilogy(timeseries.time,timeseries.vars(i,:));
end
hold off
showV([x,V],[0; timeseries.vars(:,end)],"Det. Sol. + Flds",1,nx,nx,2,1,2);

