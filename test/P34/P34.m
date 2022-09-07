clear cost
clear
%% P34 
prob.p0 = [1e-1 , 1.05 , 2.9]';% + 1e-5*randn(3,1);
prob.pbl = [0,0,0]';
prob.pbu = [100,100,10]';
prob.ibl = zeros(2,1);
prob.ib0 = 1e-1*ones(2,1);
op.tol = 1e-3;
op.rho = 1;
op.min_iter = 3;
info= SOLNP(prob,op);
cost(info.p,inf);