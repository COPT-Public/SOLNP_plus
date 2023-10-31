function f = cost(x)
   f(1) = 0;
   f(2) = 1e4*x(1)*x(2)-1;
   f(3) = exp(-x(1)) + exp(-x(2)) - 1.0001;
   f = f';
end