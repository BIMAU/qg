% Initialize QG
nx = 12;
ny = 12;
n  = nx * ny * 2;

qg = QG(nx, ny);
qg.set_par(11, 1);  % wind stress
qg.set_par(5, 45);  % Reynolds number

% load testing data
testdata = load('testdata.mat');
xr = testdata.xr(1:n);

%% Test 1: analytical Jacobian vs numerical jacobian
J = -qg.jacobian(xr, 0.0);
rhsold = qg.rhs(xr);

pert = 1e-6;
err  = zeros(n,1);
for i = 1:n
    xr(i) = xr(i) + pert;
    err(i)  = max(abs(( qg.rhs(xr) - rhsold ) / pert  - J(:,i) ));
    xr(i) = xr(i) - pert;
end

assert(norm(err,2) < 1e-4);

%% Test 2: Jacobian with time dependent contribution
dt    = 0.01;
th    = 1.0;
B     = qg.mass(n);
sig   = 1/(dt*th);
F0    = qg.rhs(xr);
pert  = 1e-6;

G = @(x) B*(x-xr)/(dt*th) + qg.rhs(x) + (1-th)/th * F0;

G0  = G(xr);

err = zeros(n,1);

J = -qg.jacobian(xr, sig);

for i = 1:n
    xr(i)  = xr(i) + pert;
    err(i) = max(abs(( G(xr) - G0 ) / pert  - J(:,i) ));
    xr(i)  = xr(i) - pert;
end

assert(norm(err,2) < 1e-4);

%% Test 3: test periodic Jacobian (analytical vs numerical)
qg = QG(nx, ny, 1);
qg.set_par(11, 1);  % wind stress
qg.set_par(5, 45);  % Reynolds number

J = -qg.jacobian(xr, 0.0);
rhsold = qg.rhs(xr);

pert = 1e-6;
err  = zeros(n,1);
for i = 1:n
    xr(i) = xr(i) + pert;
    err(i)  = max(abs(( qg.rhs(xr) - rhsold ) / pert  - J(:,i) ));
    xr(i) = xr(i) - pert;
end

assert(norm(err,2) < 1e-4);
