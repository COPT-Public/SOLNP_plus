p0 =  [10,-5]';
opt.rhoend = 1e-4;

[x,fx,oh,output] = pdfo(@fun,p0,@con,opt);
fprintf("Constraint = %e", output.constrviolation);
function f = fun(x)
    f = cost(x);
    f = f(1);
end
function [iq,eq] = con(x)
    eq = [];
    iq = cost(x);
    iq = iq(2:end);
end
% function [iq,eq] = const(x)
%     iq = [];
%     eq = [];
%     y(1) = 1.5;
%     y(2) = 2.25;
%     y(3) = 2.625;
%     for i = 1:3
%         eq = [eq;-y(i) + x(1)*(1-x(2)^i)]; 
%     end 
% end