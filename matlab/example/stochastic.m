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


m=2;
stochiter=1000;
%V=zeros(n,m);
%Y=zeros(m,stochiter);

integr_pars.ndtsub=10;
integr_pars.dt=1e-3;
integr_pars.T=10;
integr_pars.outint=0.5;

%Problem specification
%We want 
%We want a positive mass matrix since mass ought to be positive
Mass=qg.mass(n);
%Stochastic forcing
c=1e5*stochforcing(n,full(diag(Mass)));
%If the mass is positive the rhs should have a negative Jacobian
f1=@(X,b)b-qg.rhs(X); %f1(X,0) has the desired sign, we want to be able to add something to it, i.e. b.
f2=@(X)qg.jacobian(X,0.0);
%The bilinear form should have a positive sign in the rhs.
f3=@(X,Y)qg.bilin(X,Y);  
if 1
% Check on signs and consistency
    v=rand(n,1);
    ind=find(sum(spones(f2(v)),2)==1);
    v(ind)=0;
    if ( length(find(diag(Mass)<0))>0 )
       fprintf('====>>>>>> Mass matrix has wrong sign')
    end   
    if ( length(find(Mass*diag(f2(0*v))>0))>0 )
       fprintf('====>>>>>> Jacobian matrix has wrong sign')
    end 
    dif=f2(0*v)*v+2*f3(v,v)-f2(v)*v;
    if 0
    t1=f2(0*v)*v;
    t2=f3(v,v);
    t3=f2(v)*v;
    t4=f1(v,0);
    [t1([148:162]),2*t2([148:162]),-t3([148:162]),t4([148:162]) ]
    end
    norm(dif)
    w=rand(n,1);
    w(ind)=0;
    [norm(v),norm(w)]
    dif=f2(0*v)*(v+w)+2*f3(v,w)+2*f3(w,v)-f2(v)*w-f2(w)*v;
    norm(dif)
    %plot(dif)
    dif=f2(0*v)*w+f3(v,w)+f3(w,v)-f2(v)*w;
    norm(dif)
    for i=-6:3
      small=0.01^i;
      norm((f1(v+small*w,0*v)-f1(v-small*w,0*v))/(2*small) - f2(v)*w) 
    end
    pause(2)
  end
mc=size(c,2);
if (mc > m) 
  fprintf("Space stochastic forcing bigger than the stochastic space V\n");
end  
c=c(:,1);m=1;mc=1; % test with an even forcing
V=[c,randn(n,m-mc)];
showV([V],[1:m],'DO-mode',1,nx,nx,2,1,m)
Y=zeros(m,stochiter);
if 0
%test to keep x+v*[1,-1] on track
  V=mmread("v.mm");
  Y=0.1*[1,-1];
  x=mmread("x.mm");
  [norm(f1(x+V,0)), norm(f1(x-V,0)),norm(f2(x)*V)]
  %norm(f1(0*x,0))
  c=0*V;
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

