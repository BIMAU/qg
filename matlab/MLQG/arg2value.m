function [value] = arg2value(arg)
    if (ischar(arg) || isstring(arg))
        value = str2num(arg);
    else
        value = arg;
    end
end