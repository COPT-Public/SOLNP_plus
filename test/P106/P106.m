clear cost
clear
%% P106
cnstr.p0 = [5000,5000,5000,200,350,150,225,425]';% + 1e-5*randn(3,1);
cnstr.pbl = [100;1000;1000;10;10;10;10;10];
cnstr.pbu = [10000;10000;10000;1000;1000;1000;1000;1000];
cnstr.ibl = zeros(6,1);
cnstr.ib0 = 1*ones(6,1);

op.tol = 1e-4;
op.tol_con = 1e-3; 
op.qpsolver = 1;
%op.delta = 1e-1;
op.noise = 1;
op.minit = 10;
info = SOLNP(cnstr,op);
cost(info.p,inf);