function f = cost(x)
    f = -x(1)*x(2)*x(3);
    f(2) = x(1) - 4.2*sin(x(4))^2;
    f(3) =  x(2) - 4.2*sin(x(5))^2;
    f(4) =  x(3) - 4.2*sin(x(6))^2;
    f(5) = x(1) +2*x(2) + 2*x(3) - 7.2 * sin(x(7))^2;
    f = f';
     f = f.*(1+1e-4*randn(5,1));
end