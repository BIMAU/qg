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

%J = qg.jacobian(xr,0.0);
B = qg.mass(n);
assert(norm(B(:)-testdata.B(:),2) == 0);


%% Test 3: test Newton iteration

% Try to convergence from zero to 10% wind stress with Re=45:
qg.set_par(11, 0.1);  % enable  wind stress
qg.set_par(5, 45);    % 

rhs  = @ (x) qg.rhs(x);
x    = zeros(n,1);

for i = 1:20
    qg.jacob(x);
    dx = qg.solve(-qg.rhs(x));
    x = x + dx;
end

assert( norm(dx,2) < 1e-7 )

%% Test 4: perform a few backward Euler time steps

rhs = @ (x) qg.rhs(x);

x0 = zeros(n,1);

s  = 1.0/(dt*th);
B  = qg.mass(n);

F = @(x) qg.rhs(x) ;

x  = x0;
F0 = F(x);

for t = 1:3    
    for k = 1:10        

        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0; 
        qg.jacob(x, s);
        dx = qg.solve(-rhs);
        x  = x + dx;

        fprintf('%2.3e\n', norm(dx,2));
        
        if norm(dx,2) < 1e-7
            break;
        end
    end    
    fprintf('\n');
    x0 = x;
    F0 = F(x);    
end

%%% Test ...: periodic boundary conditions... todo
%
%nx = 8; ny = 8;
%n = nx * ny * 2;
%x = rand(n,1);
%qg = QG(nx, ny, 0);
%J1 = qg.jacobian(x, 0.0);
%figure(1); spy(J1);
%
%qg = QG(nx, ny, 1);
%J2 = qg.jacobian(x, 0.0);
%figure(2); spy(J2);
%figure(3); spy(J2-J1);                  