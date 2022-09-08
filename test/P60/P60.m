clear cost
clear
prob.p0 = [2,2,2]';
prob.pbl = -10 * ones(3,1);
prob.pbu = 10 * ones(3,1);
op.tol = 1e-4;
op.min_iter = 3;
op.max_iter = 30;
op.tol_con = 1e-3;
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
fprintf("f average = %e, con = %e,count = %f",f,constraint,count);
%cost(info.p,inf);