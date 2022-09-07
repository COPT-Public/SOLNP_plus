clear cost
clear subnp_dual
clear subnp_reuse
clear subnp_sr1
clear

dim = 100;
prob.p0 = ones(dim,1);% + 1e-5*randn(3,1);
prob.pbl = -10*ones(dim,1);
prob.pbu = 10*ones(dim,1);
op.tol = 1e-3;
op.tolcon = 1e-2; 
% op.maxfev = 120;
%op.delta = 1;
% op.noise = 0;
op.minit = 3;
info = SOLNP(prob,op);
cost(info.p,inf);