nx=64
Re=40
stochsize=1000
%stochsize=2
  fprintf('IS YOUR STOCHSIZE CORRECT: %4d\n', stochsize)  
%fid=fopen('timeslicesOddEvenDwsymStochForc01.txt')
% fid=fopen('timeslices.txt')
%  fid=fopen('timeslicesDt1e-3ndtsub10m4Frc1.txt')
%  fid=fopen('timeslicesDt1e-3ndtsub10m4.txt')
  fid=fopen('timeslicesForcing_dt1e-3_ndtsub10_m4.txt')
%timeslicesDt1e-3ndtsub10m4Frc10.txt')
trestart=4.5
t=trestart-1;
while t <= trestart
  [V,Y,x,t,magnitudes]=rdtimeslice(fid,stochsize);
  "read timeslice"
end
t
nx = 64;
ny = 64;
n = nx * ny * 2;
qg = QG(nx, ny);

% Set zet
zeta = 1/n;
qg.set_par(20, zeta);

%set Reynolds number
qg.set_par(5,40)
%set Full windstressfield
qg.set_par(11,1)

integr_pars.ndtsub=1;
integr_pars.dt=1e-3;
integr_pars.T=10-t;
integr_pars.outint=0.5;

Mass=qg.mass(n);
%Stochastic forcing
c=1e6*stochforcing(n,full(diag(Mass)));
%If the mass is positive the rhs should have a negative Jacobian
f1=@(X,b)b-qg.rhs(X); %f1(X,0) has the desired sign, we want to be able to add something to it, i.e. b.
f2=@(X)qg.jacobian(X,0.0);
%The bilinear form should have a positive sign in the rhs.
f3=@(X,Y)qg.bilin(X,Y);  

[x,V,Y,timings,timeseries]=DONonlin(x,V,Y,f1,f2,f3,c,Mass,integr_pars);

fclose(fid);



