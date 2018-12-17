function [y] = log_sumExp(x)
u = max(x);
y = u + log(sum(exp(x-u)));
end

