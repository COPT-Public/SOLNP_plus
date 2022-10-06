% This program use SOLNP to solve the tumor growth problem
% Input: initiate time and drug dosage
% Ouput: optimal time and drug dosage
clear
clear cost
nq = 4;
vmax = 1.1;
vcum = 65;
tmax = 200;
x = ones(1,nq)*tmax/2;
x = [x,0.5*ones(1,nq)];
prob.p0 = x';
prob.pbl = zeros(2*nq,1);
prob.pbu = [tmax*ones(nq,1);ones(nq,1)];
prob.ibl = [0;0];
prob.ibu = [vmax;vcum];
op.maxfev = 3000;
op.max_iter = 1000;
op.tol = 1e-8;
op.tol_con = 1e-8;
op.min_iter = 10;
op.delta = 1e-1;

tic
info = SOLNP_plus(prob,op);
toc
fprintf("contraint violation: %e, value = %e\n",info.constraint,info.obj);
cost(info.p,inf);

