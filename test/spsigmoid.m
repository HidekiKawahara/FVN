function output = spsigmoid(x, a, b)
output = 1.0 ./ (1 + exp(a * (x - b)));
end
