clear cost
clear
%% P34 
prob.p0 = [.1,.7,.2]';% + 1e-5*randn(3,1);
prob.pbl = [0,0,0]';
prob.ibl = 0;
prob.ib0 = 0.1;
op.tol = 1e-3;
op.tol_con = 1e-2; 
op.minit = 3;
info = SOLNP(prob,op);
cost(info.p,inf);