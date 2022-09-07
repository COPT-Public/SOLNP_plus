clear cost
clear
%% P80
prob.p0 = [-2 , 2 , 2 , -1 , -1]';% + 1e-3*randn(5,1);
prob.pbl = [-2.3,-2.3,-3.2,-3.2,-3.2]';
prob.pbu = [2.3,2.3,3.2,3.2,3.2]';
op.tol = 1e-3;
op.minit = 3;
rep = 20;
f = 0;
count = 0;
info = SOLNP(prob,op);
cost(info.p,inf);