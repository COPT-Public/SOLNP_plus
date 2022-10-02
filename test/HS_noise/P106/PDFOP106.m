% This program solve problem 106 in Hock and Schittkowski test problem
clear cost
p0 = [5000,5000,5000,200,350,150,225,425]';

lb = [100;1000;1000;10;10;10;10;10];
ub = [10000;10000;10000;1000;1000;1000;1000;1000];
opt.rhoend = 1e-3;
t = 0;
f = 0;
con = 0;
rep = 50;
count = 0;rng(1);
for i = 1:rep
    t1 = tic;
    [x,fx,oh,output] = pdfo(@fun,p0,[eye(8);-eye(8)],[ub;-lb],@const,opt);
    t = t+ toc(t1);
    con = con + constraint(x);
    f = f + fx;
    count = count + output.funcCount;
end
fprintf("time = %e,count = %f,obj = %e,con = %e\n",t/rep,count/rep,f/rep,con/rep);

function con = constraint(x)
    [iq,eq] = const(x);
    con = norm(eq)^2;
    con = con + norm(max(iq,0))^2;
    con = sqrt(con);
end
function f = fun(x)
    f = cost(x);
    f = f(1);
end
function [iq,eq] = const(x)
    iq = cost(x);
    iq = -iq(2:length(iq));
    eq = [];
end