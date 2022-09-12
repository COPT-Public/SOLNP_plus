clear cost
clear
%% P111
cnstr.p0 = -2.3*ones(10,1);% + 1e-5*randn(3,1);
cnstr.pbl = -100*ones(10,1);
cnstr.pbu = 100*ones(10,1);
op.tol = 1e-4;
op.tol_con = 1e-3; 
op.qpsolver = 1;
%op.delta = 1e-1;
op.noise = 1;
op.minit = 3;
info = SOLNP(cnstr,op);
cost(info.p,inf);