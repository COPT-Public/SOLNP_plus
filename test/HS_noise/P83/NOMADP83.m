
p0 =  ([78,33,27,27,27]' + [102,45,45,45,45]')/2;
t = 0;
con = 0;
rep = 50;
f = 0;
count = 0;
param = struct('display_degree','0','bb_output_type','OBJ PB PB PB PB PB PB','max_bb_eval',num2str(500*length(p0)),'min_poll_size','1e-4');
lb =[78,33,27,27,27]';
ub =[102,45,45,45,45]';
s = 0;
fu = @fun;
for i = 1:rep
    t1 = tic;
    [x,fval,hinf,exit_status,nfeval] = nomad(@fun,p0,lb,ub,param);
    t = t+ toc(t1);
    c = constraint(x);
    if c <= 1e-3
        s = s+1;
        f = f + fval;
        count = count + nfeval;
        con = con + c;
    end
end
fprintf("time = %e,count = %d,obj = %e,con = %e\n",t/rep,count/s,f/s,con/s)
function f = fun(x)
    f = cost(x);
    f = [f(1);-f(2:length(f));f(2:length(f))- [92,20,5]'];
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