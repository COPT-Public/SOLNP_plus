% This program solve problem 84 in Hock and Schittkowski test problem
function [f,constraint, count,t] = P84(op,rep)
prob.p0 =  [2.52,2,37.5,9.25,6.8]';% + 1e-3*randn(5,1);
prob.ibu = [294000,294000,277200]';
prob.ibl = [0,0,0]';
prob.pbu = [1000,2.4,60,9.3,7]';
prob.pbl = [0,1.2,20,9,6.5]';

% op.min_iter = 3;        
% op.tol = 1e-4;
% op.ls_time = 5;
% op.tol_con = 1e-4;rep = 1;

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