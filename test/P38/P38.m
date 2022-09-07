clear cost
clear
prob.p0 = [-3,1,-3,1]';% + 1e-4*randn(4,1);
prob.pbl = -10 * ones(4,1);
prob.pbu =  10 * ones(4,1);
op.tol = 1e-3;
op.min_iter = 3;
op.ls_time = 5; 
%op.delta = 1;
op.noise = 1;
rep = 1;
f = 0;
constraint = 0;
count = 0;
for i = 1:rep   
    info= SOLNP(prob,op);
    f = f + info.obj;
    constraint = constraint + info.constraint;
    count = count + info.count_cost;
end
f = f/rep;
count = count/rep;
constraint = constraint/rep;
fprintf("f average = %e, con = %e,count = %f\n",f,constraint,count);
%cost(info.p,inf);