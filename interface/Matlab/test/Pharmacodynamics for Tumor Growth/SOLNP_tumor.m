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
op.rescue = 1;
% op.delta = 1e-1;
t = 0;
fun.cost = @(x,nfeval)cost_para(x,nfeval,2,vmax,vcum,tmax);
rep = 1;
for i = 1:rep
    t1 = tic;
    info = SOLNP_plus(prob,op,fun);
    t = t + toc(t1);
end
fprintf("contraint violation: %e, value = %e\n",info.constraint,info.obj);
fprintf("Average running time %e\n",t/rep);
% tic
% f = cost(prob.p0);
% toc
function f = cost_para(x,nfeval,nc,vmax,vcum,tmax)
    f = zeros(nc+1,nfeval);
    for i = 1:nfeval
        f(:,i) = cost(x(:,i),1,vmax,vcum,tmax);
    end
end
