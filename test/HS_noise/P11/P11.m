% This program solve problem 11 in Hock and Schittkowski test problem
function [f,constraint, count,t,s] = P11(op,rep)
%% P11
prob.p0 = [4.9 , .1]';
prob.ib0 = ones(1,1);
prob.ibl = zeros(1,1);
% op.tol = 1e-3;
% op.tol_con = 1e-3;
% op.re_time = 0;op.ls_time = 1000;op.noise = 0;rng(1);
% op.qpsolver = 1;
% 
% %op.delta = 1e-1;
% op.min_iter = 3;
% rep = 50;rng(1);
f = 0;
constraint = 0;
count = 0;t = 0;
s = 0;rng(1);
for i = 1:rep   
    t1 = tic;
    info= SOLNP_plus(prob,op);
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
end