classdef QG < handle
    properties
        instance
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
            fprintf('QG successfully initialized\n');
        end

        function y = rhs(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_rhs(h.instance, x);
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

        function  A=jacobian(h, x, sig)
            if nargin ~= 3
                error('Two input arguments required');
            end
            [beg,jco,co]=QG_jacobian(h.instance, x, sig);
            nnz=beg(length(x)+1);
            A = crs2sp(1+double(beg),1+double(jco(1:nnz)),-co(1:nnz));
        end
       
        function  M=mass(h,n)
            if nargin ~= 2
                error('One input argument required');
            end
            M=QG_mass(h.instance,n);
            M=-spdiags(M,0,n,n);
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

        function delete(h)
            if ~isempty(h.instance)
                QG_free(h.instance);
                h.instance = [];
                fprintf('QG successfully deleted\n');
            end
        end
    end
end
