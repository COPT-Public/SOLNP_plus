clear cost
clear
%% P104
cnstr.p0 = [6,3,.4,.2,6,6,1,.5]';% + 1e-5*randn(3,1);
cnstr.pbl = .1*ones(8,1);
cnstr.pbu = 10*ones(8,1);
cnstr.ibl = [zeros(4,1);1];
cnstr.ibu = [inf*ones(4,1);4.2];
cnstr.ib0 = 3*ones(5,1);
op.tol = 1e-4;
op.tol_con = 1e-4; 
op.qpsolver = 1;
%op.delta = 1e-1;
op.noise = 1;
op.minit = 3;
info = SOLNP(cnstr,op);
cost(info.p,inf);