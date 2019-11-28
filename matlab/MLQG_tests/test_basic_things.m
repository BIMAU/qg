% Initialize QG
nx = 32;
ny = 32;
n  = nx * ny * 2;
qg = QG(nx, ny);
qg.set_par(11, 0.0);  % wind stress
qg.set_par(5, 45);    % Reynolds number

%% Test 1: rhs

x = zeros(n, 1);
F = qg.rhs(x);
assert(numel(F) == n);
assert(norm(F,2) == 0);


%% Test 2: test matrix and rhs against data

testdata = load('testdata.mat');
qg.set_par(11, 0.1);  % enable some wind stress
x =(-3+mod(1:n,7))/7; % some nonzero entry
y = qg.apply(x);
assert(norm(y(:)-testdata.y(:),2) == 0);

F = qg.rhs(x);
assert(norm(F(:)-testdata.F(:),2) == 0);

xr = testdata.xr;
yr = qg.apply(xr);
assert(norm(yr(:)-testdata.yr(:),2) == 0);

J = qg.jacobian(xr,0.0);