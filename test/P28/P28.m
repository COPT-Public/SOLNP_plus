clear cost
clear
%% P28
prob.p0 = [-4,1,1]'; %+ 1e-4*randn(3,1);
op.tol = 1e-3;
op.min_iter = 3;
op.tol_con = 1e-2;
% op.noise = 0;
rep = 1;
f = 0;
t = 0;
constraint = 0;
count = 0;
for i = 1:rep   
    t1 = tic;
    info= SOLNP(prob,op);
    t = t + toc(t1);
    f = f + info.jh(length(info.jh));
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