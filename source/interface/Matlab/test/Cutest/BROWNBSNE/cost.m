function f = cost(x)
   f(1) = 0;
   f(2) = x(1) - 1e6;
   f(3) = x(2)- 2*1e-6;
   f(4) = x(1)*x(2) - 2;
   f = f';
end