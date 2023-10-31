% This program solve problem 11 in Hock and Schittkowski test problem
clear cost
p0 = [-3,1,-3,1]';
opt.rhoend = 1e-3;
opt.solver = "cobyla";
t = 0;
con = 0;
rep = 50;
f = 0;
count = 0;
rng(1);
for i = 1:rep
    t1 = tic;
    [x,fx,oh,output] = pdfo(@fun,p0,[eye(4);-eye(4)],[10 * ones(4,1);10 * ones(4,1)],@const,opt);
    t = t+ toc(t1);
    f = f +fx;
    count = count + output.funcCount;
    con = con + constraint(x);
end
fprintf("time = %e,count = %f,obj = %e,con = %e\n",t/rep,count/rep,f/rep,con/rep)

function f = fun(x)
    f = cost(x);
    f = f(1);
end
function con = constraint(x)
    [iq,eq] = const(x);
    con = norm(eq)^2;
    con = con + norm(max(iq,0))^2;
    con = sqrt(con);
end
function [iq,eq] = const(x)
    eq = cost(x);
    eq = eq(2:length(eq));
    iq = [];
end