% This program solve problem 11 in Hock and Schittkowski test problem
% function [f,constraint, count,t] = P11(op,rep)
%% P11
clear
prob.p0 = [10,-5]';
prob.ib0 = -ones(3,1);
prob.ibu = zeros(3,1);

op.ls_way = 1;
op.min_iter = 5;
op.noise = 1;
op.tol = 1e-4;
op.ls_time = 5;
op.tol_con = 1e-4;rep = 1;
op.tol_restart = 1e0;
% op.re_time = 10;
fun.cost = @(x,nfeval)cost_para(x,nfeval);
% fun.grad = @grad;
% fun.hess = @hess;

f = 0;
constraint = 0;
count = 0;t = 0;
s = 0;rng(1);
for i = 1:rep   
    t1 = tic;
    info= SOLNP(prob,op,fun);
    t = t + toc(t1);
        s = s +1;
        f = f +  info.obj;
        constraint = constraint + info.constraint;
        count = count + info.count_cost;
    clear cost
end
t = t/rep;
f = f/s;
count = count/s;
constraint = constraint/s;
fprintf("f average = %e, con = %e,count = %f,time = %e\n",f,constraint,count,t);
%cost(info.p,inf);
% end

function f = cost_para(x,nfeval)
    f = cost(x(:,1));
    for i = 2:nfeval
        f(:,i) = cost(x(:,i));
    end
end
