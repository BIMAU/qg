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

%% Test 2: test matrices and rhs against data
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
B = qg.mass(n);
assert(norm(B(:)-testdata.B(:),2) == 0);

%% Test 3: test Newton iteration
% Try to converge from zero to 10% wind stress with Re=45:
qg.set_par(11, 0.01); % enable wind stress
qg.set_par(5, 45);    %

x   = zeros(n,1);
fprintf('\nNewton\n')
for i = 1:20
    qg.jacob(x);
    qg.compute_precon();
    dx = qg.solve(-qg.rhs(x));
    x  = x + dx;
    fprintf('|dx| = %e\n', norm(dx,2))
    if norm(dx,2) < 1e-7
        break;
    end
end
assert( i < 20 )
assert( norm(dx,2) < 1e-7 )

%% Test 4: perform a few backward Euler time steps
rhs = @ (x) qg.rhs(x);

x0 = zeros(n,1);

dt = 0.01;
th = 1;
s  = 1.0/(dt*th);
B  = qg.mass(n);

F = @(x) qg.rhs(x) ;

x  = x0;
F0 = F(x);

fprintf('\nTimestepping\n')
for t = 1:3
    for k = 1:10
        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
        qg.jacob(x, s);
        qg.compute_precon();
        dx = qg.solve(-rhs);
        x  = x + dx;
        fprintf('t = %d, k = %d, |dx| = %e\n', t, k, norm(dx,2))
        if norm(dx,2) < 1e-7
            break;
        end
    end
    assert(norm(dx,2) < 1e-7);
    x0 = x;
    F0 = F(x);
end