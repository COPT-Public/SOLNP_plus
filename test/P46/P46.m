clear cost
clear
%% P46
prob.p0 = [.5*sqrt(2),1.75,.5,2,2]';% + 1e-3*randn(5,1);
op.tol = 1e-3;
op.ls_time = 3;
op.min_iter = 4;
op.tolcon = 1e-2;
rep = 1;
f = 0;
constraint = 0;
count = 0;t = 0;
for i = 1:rep   
    t1 = tic;
    info= SOLNP(prob,op);
    t = t + toc(t1);
    f = f + + info.obj;
    constraint = constraint + info.constraint;
    count = count + info.count_cost;
    clear cost
end
t = t/rep;
f = f/rep;
count = count/rep;
constraint = constraint/rep;
fprintf("f average = %e, con = %e,count = %f,time = %e\n",f,constraint,count,t);
%cost(info.p,inf);