clear cost
clear
prob.p0 = [0,0]';
op.tol = 1e-3;
op.tolcon = 1e-2;
op.min_iter = 3;
rep = 1;
f = 0;
constraint = 0;
count = 0;
t = 0;
for i = 1:rep   
    t1 = tic;
    info= SOLNP_qp(prob,op);
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
fprintf("f average = %e, con = %e,count = %f,time = %f",f,constraint,count,t);
%cost(info.p,inf);