% This program solve problem 84 in Hock and Schittkowski test problem

p0 = [2.52,2,37.5,9.25,6.8];
opt.rhoend = 1e-4;
t = 0;
f = 0;
con = 0;
rep = 1;
count = 0;
s = 0;rng(1);
for i = 1:rep
    t1 = tic;
    [x,fx,oh,output] = pdfo(@fun,p0,[-eye(5);eye(5)],[ -[0,1.2,20,9,6.5]';[1000,2.4,60,9.3,7]'],@const,opt);
    td = toc(t1);
    c = constraint(x);
    if c <=1e-3
        t = t+ td;
        s = s+1;
        con = con + c;
        f = f + fx;
        count = count + output.funcCount;
    end
end

fprintf("time = %e,count = %f,obj = %e,con = %e\n",t/s,count/s,f/s,con/s);
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
    iq = [-iq;iq -  [294000,294000,277200]'];
    eq = [];
end