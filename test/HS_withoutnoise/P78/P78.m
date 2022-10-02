% This program solve problem 78 in Hock and Schittkowski test problem
function [f,constraint, count,t] = P78(op,rep)
%% P78
prob.p0 =[-2,1.5,2,-1,-1]';%+1e-4*randn(5,1);
% op.tol = 1e-3;
% op.tol_con = 1e-3;
% op.min_iter = 3;
% % op.ls_time = 5;
% % op.noise = 0;
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