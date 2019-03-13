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