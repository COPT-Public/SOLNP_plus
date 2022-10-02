% This program solve problem 93 in Hock and Schittkowski test problem
function [f,constraint, count,t] = P93(op,rep)
%% P93
prob.p0 = [5.54,4.4,12.02,11.82,.702,.852]';% + 1e-5*randn(3,1);
prob.pbl = zeros(6,1);
prob.ibl = [0;0];
prob.ib0 = [1;1];

% op.min_iter = 3;        
% op.tol = 1e-4;
% op.ls_time = 5;
% op.tol_con = 1e-4;rep = 1;
% 

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
fprintf("f average = %e, con = %e,count = %f,time = %e\n, s= %d",f,constraint,count,t,s);
%cost(info.p,inf);
end