function f = cost(x)
   f(1) = -x(1) + 10*(x(1)^2+x(2)^2-1);
   f(2) = x(1)^2+x(2)^2-1;
   f = f';
end