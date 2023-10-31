% This program solve problem 83 in Hock and Schittkowski test problem
p0 = ([78,33,27,27,27]' + [102,45,45,45,45]')/2;
opt.rhoend = 1e-3;
t = 0;
f = 0;
con = 0;
rep = 50;
count = 0;
s = 0;rng(1);
for i = 1:rep
    t1 = tic;
    [x,fx,oh,output] = pdfo(@fun,p0,[eye(5);-eye(5)],[[102,45,45,45,45]';-[78,33,27,27,27]'],@const,opt);
    t = t+ toc(t1);
    c = constraint(x);
    if c <=1e-3
        s = s+1;
        con = con + c;
        f = f + fx;
        count = count + output.funcCount;
    end
   
end

fprintf("time = %e,count = %f,obj = %e,con = %e\n",t/rep,count/s,f/s,con/s);
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
    iq = cost(x);
    iq = iq(2:length(iq));
    iq = [-iq;iq -  [92,20,5]'];
    eq = [];
end