clear cost
clear
%% P93
cnstr.p0 = [5.54,4.4,12.02,11.82,.702,.852]';% + 1e-5*randn(3,1);
cnstr.pbl = zeros(6,1);
cnstr.ibl = [0;0];
cnstr.ib0 = [1e-2;1e-2];
op.tol = 1e-3;
op.tol_con = 1e-2; 
op.qpsolver = 1;
%op.maxfev = 120;
%op.delta = 1;
op.noise = 1;
op.minit = 3;
info = SOLNP(cnstr,op);
cost(info.p,inf);