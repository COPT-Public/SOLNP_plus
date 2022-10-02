% This program solve problem 38 in Hock and Schittkowski test problem
function [f,constraint, count,t] = P38(op,rep)
prob.p0 = [-3,1,-3,1]';
prob.pbl = -10 * ones(4,1);
prob.pbu =  10 * ones(4,1);
% op.tol = 1e-3;
% % op.ls_time = 5;
% op.tol_con = 1e-3;
% %op.delta = 1;
% op.min_iter = 3;
% op.re_time = 0;op.ls_time = 1000;op.noise = 0;rng(1);
% rep = 50;
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
fprintf("f average = %e, con = %e,count = %f,time = %e\n",f,constraint,count,t);
%cost(info.p,inf);
end