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
  
x = zeros(n, 1);
tol=1e-10;


m=3;
stochiter=1000;
V=zeros(n,m);
Y=zeros(m,stochiter);
Mass=qg.mass(n);
%vsm(Mass)

integr_pars.ndtsub=1;
integr_pars.dt=1e-2;
integr_pars.T=10;
integr_pars.outint=0.5;

c=1e6*stochforcing(n,full(diag(Mass)));

%c=1e-5*ones(n,1);
f1=@(X,b)b-qg.rhs(X);
f2=@(X)qg.jacobian(X,0.0);
f3=@(X,Y)-qg.bilin(X,Y); 
 if 1
    v=rand(n,1);
    ind=find(sum(spones(f2(v)),2)==1);
    v(ind)=0;
    %vsm(f2(0*v));
    dif=f2(0*v)*v-2*f3(v,v)-f2(v)*v;
    if 0
    t1=f2(0*v)*v;
    t2=f3(v,v);
    t3=f2(v)*v;
    t4=f1(v,0);
    [t1([148:162]),-2*t2([148:162]),-t3([148:162]),t4([148:162]) ]
    end
    norm(dif)
    w=rand(n,1);
    w(ind)=0;
    [norm(v),norm(w)]
    dif=f2(0*v)*(v+w)-2*f3(v,w)-2*f3(w,v)-f2(v)*w-f2(w)*v;
    norm(dif)
    %plot(dif)
    dif=f2(0*v)*w-f3(v,w)-f3(w,v)-f2(v)*w;
    norm(dif)
    for i=-6:3
      small=0.01^i;
      norm((f1(v+small*w,0*v)-f1(v-small*w,0*v))/(2*small) - f2(v)*w) 
    end
    pause(2)
  end 
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
showV([x,V],[0; timeseries.vars(:,end)],"Det. Sol. + Flds",1,nx,nx,2,1,4);

