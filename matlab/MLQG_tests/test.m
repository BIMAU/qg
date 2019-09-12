nx = 128;
ny = 128;
n  = nx * ny * 2;
idx1 = 1:2:n;
idx2 = 2:2:n;

qg = QG(nx, ny);

% Set zeta
zeta = 1/n;
qg.set_par(20, zeta);

qg.set_par(11, 0.0101);

x = zeros(n, 1);
x = newton(qg, x, 1e-10);

xhopf    = load('fort.12');
xhopflin = load('fort.12.lin');
Fhopf    = load('fort.13');

for i = 1:21
    fprintf('par %d: %1.8f\n', i, qg.get_par(i));
end
F = qg.rhs(xhopf);

figure(1)
contourf(reshape(x(idx2),nx,ny)');
colorbar

figure(2)
contourf(reshape(xhopf(idx2),nx,ny)');
colorbar

J = qg.jacobian(xhopf,0.0);
J = J([idx1,idx2],[idx1,idx2]);

% Load Hopf matrix
beg  = importdata('C.beg');
co   = importdata('C.co');
jco  = importdata('C.jco');

Jhopf = crs2sp(beg,jco,co);
Jhopf = Jhopf([idx1,idx2],[idx1,idx2]);