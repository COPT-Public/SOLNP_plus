
% This program solve 
% min_x ||x||^2
% subject to -10<=x<=10
clear cost
clear

dim = 100;
prob.p0 = ones(dim,1);% + 1e-5*randn(3,1);
prob.pbl = -10*ones(dim,1);
prob.pbu = 10*ones(dim,1);
op.tol = 1e-4;
op.tol_con = 1e-4; 
op.min_iter = 10;

info = SOLNP(prob,op);
cost(info.p,inf);