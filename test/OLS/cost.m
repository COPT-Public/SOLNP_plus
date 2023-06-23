
% function f = cost(x)
%     con = 40;
%     f = cost_in(x(con+1:end));
%     f(1) = sum(x(1:con).^2);
%     f(2:end) = x(1:con) - f(2:end);
% end

function f = cost(x)
    rng(223)
    a = [2,3,-5];
    num = 2;
    t = linspace(-3,3,num);
    y = a(1) +a(2)*t+ a(3)*t.^2; %+ 1e-4*randn(num,1);
%     y = a(1) + a(2)*t + exp(a(3)*t);
    f = 0;
    for i = 1:num
%         f = [f;(y(i)-(x(1)+x(2)*t(i)+ exp(x(3)*t(i))))^2];
%         f = f + (y(i)-(x(1)+x(2)*t(i)+ x(3)*t(i)^2))^2;
         f = [f;(y(i)-(x(1)+x(2)*t(i)+ x(3)*t(i)^2))^2];
    end
end
