
p0 = [-2.6 , 2 ,2]';
t = 0;
con = 0;
rep = 50;
f = 0;
count = 0;
param = struct('display_degree','0','bb_output_type','OBJ PB PB','max_bb_eval',num2str(500*length(p0)),'min_poll_size','1e-3');
lb = -inf*ones(3,1);
ub = -lb;
rng(1);
for i = 1:rep
    t1 = tic;
    [x,fval,hinf,exit_status,nfeval] = nomad(@fun,p0,lb,ub,param);
    t = t+ toc(t1);
    f = f + fval;
    count = count + nfeval;
    con = con + constraint(x);
end
fprintf("time = %e,count = %d,obj = %e,con = %e\n",t/rep,count/rep,f/rep,con/rep)
function f = fun(x)
    f = cost(x);
    f = [f;-f(2:length(f))];
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