% Initialize QG
nx = 32;
ny = 32;
n  = nx * ny * 2;
qg = QG(nx, ny);
qg.set_par(11, 0.0);  % Wind stress
qg.set_par(5, 45);    % Re number

%% Test 1: RHS
x = zeros(n, 1);
F = qg.rhs(x);
assert(numel(F) == n);

%% Test 2: matrix action
x = rand(n, 1);
