classdef QG < handle
    properties
        instance

        % grid points in x-direction
        nx

        % grid points in y-direction
        ny

        % number of unknowns
        nun = 2

        % mass matrix
        B

        % theta in theta-method time stepping #TODO
        theta = 1.0

        % maximum Newton iterations
        Nkmx = 10

        % Newton tolerance
        Ntol = 1e-3
    end
    methods
        function h = QG(nx, ny, perio)
            if ((nargin < 2) || (nargin > 3))
                error('Wrong number of input arguments');
            end
            if (nargin ~= 3)
                perio = 0;
            end
            h.instance = QG_init(nx, ny, perio);
            h.nx       = nx;
            h.ny       = ny;
            h.B        = h.mass(h.nx * h.ny * h.nun);
            fprintf('QG successfully initialized\n');
        end

        function y = rhs(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_rhs(h.instance, x);
        end

        function [u,v] = compute_uv(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            [u,v] = QG_compute_uv(h.instance, x(:));
        end
        
        function [gradx, grady] = gradient(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            [gradx, grady] = QG_gradient(h.instance, x(:));
        end

        function [x, y] = grid(h)
            if nargin ~= 1
                error('No input argument required');
            end
            [x, y] = QG_grid(h.instance);
        end

        function Z = bilin(h, V,W)
            if nargin ~= 3
                error('Two input arguments required');
            end
            mV=size(V,2);mW=size(W,2);
            for i=1:mV
                for j=1:mW
                    Z(:,(i-1)*mW+j) = QG_bilin(h.instance, V(:,i),W(:,j));
                end
            end
        end

        function jacob(h, x, sig)
            if ((nargin < 2) || (nargin > 3))
                error('Wrong number of input arguments');
            end
            if nargin < 3
                sig = 0.0;
            end
            QG_jacob(h.instance, x, sig);
        end

        function  A = jacobian(h, x, sig)
            if nargin ~= 3
                error('Two input arguments required');
            end
            [beg,jco,co]=QG_jacobian(h.instance, x, sig);
            nnz=beg(length(x)+1);
            A = crs2sp(1+double(beg),1+double(jco(1:nnz)),-co(1:nnz));
        end

        function  M = mass(h,n)
            if nargin ~= 2
                error('One input argument required');
            end
            M = QG_mass(h.instance,n);
            M = -spdiags(M,0,n,n);
        end

        function y = apply(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_apply(h.instance, x);
        end

        function y = compute_precon(h)
            if nargin ~= 1
                error('No input arguments required');
            end
            QG_compute_precon(h.instance);
        end

        function y = solve(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_solve(h.instance, x);
        end

        function set_par(h, par, val)
            if nargin ~= 3
                error('Two input arguments required');
            end
            QG_set_par(h.instance, par, val);
        end

        function val = get_par(h, par)
            if nargin ~= 2
                error('One input argument required');
            end
            val = QG_get_par(h.instance, par);
        end

        function [x, k] = step(h, x, dt)
        % time step

            xm   = x;
            rhsm = h.rhs(xm);
            s    = 1.0/(dt*h.theta);

            % Newton
            for k = 1:h.Nkmx
                rhs = h.B*(x-xm)*s + h.rhs(x) ...
                      + (1-h.theta)/h.theta * rhsm;
                jac = h.jacobian(x, s);
                dx  = jac \ rhs;
                x   = x + dx;

                if norm(dx,2) < h.Ntol
                    break;
                end
            end

            if k == h.Nkmx
                ME = MException('QG:convergenceError', ...
                                'no convergence in Newton iteration');
                throw(ME);
            end
        end

        function delete(h)
            if ~isempty(h.instance)
                QG_free(h.instance);
                h.instance = [];
                fprintf('QG successfully deleted\n');
            end
        end
    end
end
