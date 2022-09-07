clear cost
clear
%% P100
cnstr.p0 = [1,2,0,4,0,1,1]';% + 1e-5*randn(3,1);
cnstr.ibl = zeros(4,1);
cnstr.ib0 = 1*ones(4,1);
op.tol = 1e-4;
op.tol_con = 1e-2; 
op.qpsolver = 1;
%op.delta = 1e-1;
% op.maxfev = 200;
op.noise = 0;
op.minit = 3;
info = SOLNP(cnstr,op);
cost(info.p,inf);