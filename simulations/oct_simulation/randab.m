function out = randab(a, b)
% This is a function to generate random number between two enterval
m = min(a, b);
M = max(a, b);
out = m + rand(1) * (M-m);
end
