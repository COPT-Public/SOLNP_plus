clear cost
clear
prob.p0 = [-1.2,1]';
op.tol = 1e-3;
%op.tolcon = 1e-2;
op.min_iter = 3;
op.tol_con = 1e-2;
%op.maxfev = 47;
% op.noise = 0;
rep = 20;
f = 0;
constraint = 0;
count = 0;
t = 0;
for i = 1:rep   
    t1 = tic;
    info= SOLNP(prob,op);
    t = t + toc(t1);
    f = f + info.obj;
    constraint = constraint + info.constraint;
    count = count + info.count_cost;
    clear cost
end
t = t/rep;
f = f/rep;
count = count/rep;
constraint = constraint/rep;
fprintf("f average = %e, con = %e,count = %f,time = %f\n",f,constraint,count,t);
%cost(info.p,inf);