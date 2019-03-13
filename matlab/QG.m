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

        function delete(h)
            if ~isempty(h.instance)
                QG_free(h.instance);
                h.instance = [];
                fprintf('QG successfully deleted\n');
            end
        end
    end
end