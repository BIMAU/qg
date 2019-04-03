classdef QG < handle
    properties
        instance
    end
    methods
        function h = QG(A, params)
            if nargin ~= 2
                error('Two input arguments required');
            end
            h.instance = QG_init(A, params);
            fprintf('QG successfully initialized\n');
        end

        function y = rhs(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_rhs(h.instance, x);
        end
        
	function z = bilin(h, x,y)
            if nargin ~= 3
                error('Two input arguments required');
            end
            z = QG_bilin(h.instance, x,y);
        end

        function jacob(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            QG_jacob(h.instance, x);
        end

        function  A=jacobian(h, x, sig)
            if nargin ~= 3
                error('Two input arguments required');
            end
            [beg,jco,co]=QG_jacobian(h.instance, x, sig);
            nnz=beg(length(x)+1);
           % size(beg)
	   % A = csr2sp (beg+1,jco+1,co);
           A = csr2sp (1+double(beg),1+double(jco(1:nnz)),-co(1:nnz));
        end
       
        function  M=mass(h,n)
            if nargin ~= 2
                error('One input argument required');
            end
            M=QG_mass(h.instance,n);
            M=-spdiags(M,0,n,n)
        end

        function y = apply(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            y = QG_apply(h.instance, x);
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
