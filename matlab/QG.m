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

        function jacob(h, x)
            if nargin ~= 2
                error('One input argument required');
            end
            QG_jacob(h.instance, x);
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