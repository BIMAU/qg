function [] = print0(varargin)
    global pid
    if pid == 0
        fprintf('pid0: ');
        fprintf(varargin{:});
    end
end