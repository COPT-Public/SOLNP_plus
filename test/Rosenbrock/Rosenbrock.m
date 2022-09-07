clear cost
clear subnp_dual
clear subnp_reuse
clear subnp_sr1
clear

prob.p0 = zeros(2,1);% + 1e-5*randn(3,1);
prob.pbl = -1*ones(2,1);
prob.pbu = 2*ones(2,1);
op.tol = 1e-4;
op.tolcon = 1e-2; 
% op.maxfev = 120;
%op.delta = 1;
% op.noise = 0;
op.minit = 10;
info = SOLNP(prob,op);
cost(info.p,inf);