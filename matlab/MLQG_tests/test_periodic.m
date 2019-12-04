% Initialize QG
nx = 12;
ny = 12;
n  = nx * ny * 2;
qg = QG(nx, ny, 1);

qg.set_par(11, 0.01);  % wind stress
qg.set_par(5, 45);  % Reynolds number

% load testing data
testdata = load('testdata.mat');
xr = testdata.xr(1:n);

%% Test 1: Newton

qg.set_par(11, 0.001); % enable  wind stress
qg.set_par(5, 45);     % 

rhs = @ (x) qg.rhs(x);
x   = zeros(n,1);

for i = 1:2
    qg.jacob(x);
    dx = qg.solve(-qg.rhs(x));
    x  = x + dx;
    fprintf('%f\n', norm(dx,2));
end

%assert( norm(dx,2) < 1e-7 )

%% Test 2: analytical jacobian vs numerical jacobian
J = -qg.jacobian(xr, 0.0);
rhsold = qg.rhs(xr);

pert = 1e-6;
err  = zeros(n,1);
for i = 1:n
    xr(i) = xr(i) + pert;
    err(i)  = max(abs(( qg.rhs(xr) - rhsold ) / pert  - J(:,i) ));
    xr(i) = xr(i) - pert;
end

%assert(norm(err,2) < 1e-4);

plot(err)
