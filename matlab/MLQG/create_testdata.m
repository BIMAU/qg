% Initialize QG
nx = 32;
ny = 32;
n  = nx * ny * 2;
qg = QG(nx, ny);

qg.set_par(5, 45);    % Reynolds number
qg.set_par(11, 0.1);  % enable some wind stress
x =(-3+mod(1:n,7))/7; % some nonzero entry

y  = qg.apply(x);
F  = qg.rhs(x);
xr = testdata.xr;
yr = qg.apply(xr);

J = qg.jacobian(xr,0.0);
B = qg.mass(n);
