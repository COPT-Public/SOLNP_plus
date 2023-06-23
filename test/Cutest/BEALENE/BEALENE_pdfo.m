p0 =  [1,1]';
opt.rhoend = 1e-4;

[x,fx,oh,output] = pdfo(@fun,p0,opt);

function f = fun(x)
    f = 0;
    y(1) = 1.5;
    y(2) = 2.25;
    y(3) = 2.625;
    for i = 1:3
        f = f + (-y(i) + x(1)*(1-x(2)^i))^2; 
    end 
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