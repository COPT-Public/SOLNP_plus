a = asin(sqrt(1/4.2));
b = asin(sqrt(5/7.2));
p0 = [ones(1,3),a,a,a,b]';
opt.rhoend = 1e-3;
t = 0;f=0;count = 0;
rep = 50;s = 0;con = 0;
rng(1);
for i = 1:rep
    t1 = tic;
    [x,fx,oh,output] = pdfo(@fun,p0,@const,opt);
    c = constraint(x);
    if c <= 1e-3
        con = con + c;
        count = count + output.funcCount;
        s = s+1;
        f = f+fx;
    end
    t = t+ toc(t1);
end

fprintf("time = %e,count = %2f,obj = %e,con = %e\n",t/rep,count/s,f/s,con/s)
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