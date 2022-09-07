clear cost
clear
%% P2  
prob.p0 = [-2.6 , 2 ,2]';% + 1e-5*randn(3,1);
op.tol = 1e-2;
op.tol_con= 1e-2;
op.rho = 1;
op.noise = 1;
op.qpsolver = 1;
op.min_iter = 3;
info= SOLNP(prob,op);
cost(info.p,inf);