% Initialize QG
nx = 16;
ny = 16;
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

norm(err,2)

%assert(norm(err,2) < 1e-4);

%vsm(J)
%vsm(Jn-J)

%% Test 2: Newton convergence
qg.set_par(11, 0.01); % enable  wind stress
qg.set_par(5, 25);     % Reynolds number

rhs = @ (x) qg.rhs(x);
x   = zeros(n,1);

for i = 1:10
    qg.jacob(x);
    dx = qg.solve(-qg.rhs(x));
    x  = x + dx;
    fprintf('%2.13f\n', norm(qg.rhs(x),2));
end
 
%assert( norm(dx,2) < 1e-4 )