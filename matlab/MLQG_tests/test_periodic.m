% Initialize QG
nx = 32;
ny = 32;
n  = nx * ny * 2;
qg = QG(nx, ny, 1);

qg.set_par(11, 0.1);  % wind stress
qg.set_par(5, 45);    % Reynolds number

% load testing data
testdata = load('testdata.mat');
xr = testdata.xr(1:n);

%% Test 1: analytical jacobian vs numerical jacobian

J  = -qg.jacobian(xr, 0.0);

F0 = qg.rhs(xr);

pert = 1e-6;
err  = zeros(n,1);
D    = zeros(n,n);
Jn   = zeros(n,n);

v = zeros(n,1);
for i = 1:n
    v(i) = pert;
    Jn(:,i) = (qg.rhs(xr+v) - qg.rhs(xr-v)) / (2*pert);
    v(i) = 0;
    
    D(:,i)  = abs(Jn(:,i)  - J(:,i));
    err(i)  = max(D(:,i));
end

% --> check dit in relatieve zin...
Jn = sparse(Jn);
D(abs(D) < 1) = 0.0;
D = sparse(D);

assert(max(err) < 1e-4);

%vsm(J)
%vsm(Jn-J)

%% Test 2: perform a few backward Euler time steps
qg.set_par(11, 0.1); % enable  wind stress
qg.set_par(5, 45); % Reynolds number

rhs = @ (x) qg.rhs(x);

x0 = zeros(n,1);

dt = 0.001;
th = 1;
s  = 1.0/(dt*th);
B  = qg.mass(n);

F = @(x) qg.rhs(x) ;

x  = x0;
F0 = F(x);

kDes = 3.5;
for t = 1:200
    fprintf('dt = %2.2e, Newton: \n', dt);
    plotQG(nx,ny,2,x)
    drawnow
        
    for k = 1:10        
        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0; 
        qg.jacob(x, s);
        dx = qg.solve(-rhs);
        x  = x + dx;
        fprintf('||dx|| = %2.5e\n', norm(dx));
        if norm(dx,2) < 1e-7
            break;
        end
    end
    
    dt = kDes / k * dt;
    s  = 1.0/(dt*th);
    
    assert(norm(dx,2) < 1e-7);
    
    x0 = x;
    F0 = F(x);    
end