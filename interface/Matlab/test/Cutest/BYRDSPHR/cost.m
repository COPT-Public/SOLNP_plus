function f = cost(x)
   f(1) = -x(1)-x(2)-x(3);
   f(2) = sum(x.^2)-9;
   f(3) = (x(1)-1)^2+x(2)^2+x(3)^2-9;
   f = f';
end