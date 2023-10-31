function f = cost(x)
   f(1) = -x(1);
   f(2) = x(2)- x(1)^3-x(3)^2;
   f(3) = x(1)^2 - x(2) - x(4)^2;
   f = f';
end