function [exp_x ,max_x] = NormalizedExp(x)
max_x = max(x);
x = x-max_x;
exp_x = exp(x);