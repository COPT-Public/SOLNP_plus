clear
% This program solve problem 106 in Hock and Schittkowski test problem
%% P106
% function [f,constraint, count,t] = P106(op,rep)
prob.p0 = [5000,5000,5000,200,350,150,225,425]';% + 1e-5*randn(3,1);
prob.pbl = [100;1000;1000;10;10;10;10;10];
prob.pbu = [10000;10000;10000;1000;1000;1000;1000;1000];
prob.ibl = zeros(6,1);
prob.ib0 = 1*ones(6,1);

op.ls_way = 1;
op.min_iter = 10;        
op.tol = 1e-4;
op.ls_time = 10;
op.tol_con = 1e-4;rep = 1;
rng(1);
fun.cost = @(x,nfeval)cost_para(x,nfeval,6);

f = 0;
constraint = 0;
count = 0;t = 0;
s = 0;rng(1);
for i = 1:rep   
    t1 = tic;
    info= SOLNP(prob,op,fun);
    t = t + toc(t1);
    if info.constraint<= 1e-3 
        s = s +1;
        f = f +  info.obj;
        constraint = constraint + info.constraint;
        count = count + info.count_cost;
    end
    clear cost
end
t = t/rep;
f = f/s;
count = count/s;
constraint = constraint/s;
fprintf("f average = %e, con = %e,count = %f,time = %e, success = %d\n",f,constraint,count,t,s);
%cost(info.p,inf);
% end

function f = cost_para(x,nfeval,nc)
    f = zeros(nc+1,nfeval);
    for i = 1:nfeval
        f(:,i) = cost(x(:,i));
    end
end