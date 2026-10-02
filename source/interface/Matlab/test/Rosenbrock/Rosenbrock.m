
% This program minimize the Rosenbrock function in a box 
clear
clear cost

prob.p0 = [-1.2;1];%zeros(2,1);% + 1e-5*randn(3,1);
% prob.pbl = -1*ones(2,1);
% prob.pbu = 2*ones(2,1);
op.tol = 1e-4;
op.tol_con = 1e-4; 
op.min_iter = 5;
fun.cost = @(x,nfevel)cost_para(x,nfevel,0);
info = SOLNP_plus(prob,op,fun);
cost(info.p,inf);

function f = cost_para(x,nfeval,nc)
    f = zeros(nc+1,nfeval);
    for i = 1:nfeval
        f(:,i) = cost(x(:,i));
    end
end