% This program solve 
% min_x ||x||^2
% subject to -10<=x<=10

clear cost
clear
rng(13)
% prob.pbl = -10*ones(dim,1);
% prob.pbu = 10*ones(dim,1);
dim = 100;
a = abs(1+5e-1*randn(dim));
fun_l2 = @(x) norm(a*x)^2;

op.tol = 1e-6;
op.tol_con = 1e-4; 
op.noise = 1;
op.delta = 1e-5;
op.min_iter = 1;
op.max_iter = 10000;
% op.bfgs = 0;
op.max_fev = op.max_iter *3;
op.ls_time = 10;
op.re_time = op.max_iter;
% op.scale = 0;
op.tol_restart = 1e-4;

% fun.grad = @x;
fun.cost = @(x,nfeval)cost_para(x,nfeval,0,fun_l2);
num = 1;
% dim = floor(logspace(2,2.5,num));
fd_result = zeros(num,1);
rs_result = zeros(num,1);
prob.p0 = ones(dim,1);
% for i = 1:num
% %     prob.p0 = ones(dim(i),1);
%     info_full = SOLNP(prob,op,fun);
% 
% %     fd_result(i) =  fun(info.p,inf);
%     clear cost
% end

% loglog(dim,fd_result);
% hold on

op.tol = 1e-6;
op.bfgs = 0;
op.scale = 0;
op.noise = 0;
op.drsom = 1;
op.delta = 1e-4;
op.rs = 0;count = 0;obj = 0;
% fun.grad = @x;
rep = 1;
for i = 1:num
    op.batchsize = 20;%floor(dim/4);
%     prob.p0 = ones(dim(i),1);
    count = 0;
    for j = 1:rep
        info_rs = SOLNP(prob,op,fun,0,zeros(dim(i)));
        count = count + info_rs.count_cost;
        obj = obj + info_rs.obj;
%         temp = cost(info.p,inf);
%         if temp < 1e7
%             rs_result(i) = rs_result(i) + temp;
%             count = count + 1;
%         end
%             clear cost
    end
    count = count/rep;
    obj = obj/rep;
%     if count > 0
%         rs_result(i) = rs_result(i)/count;
%     else
%         rs_result(i) = inf;
%     end
end
% loglog(dim,rs_result);
% hold on
% legend("Finite difference","Random Sampling");
% ylabel("Function Evaluation");
% xlabel("Dimension of the Problem")

function f = cost_para(x,nfeval,nc,fun)
    f = zeros(nc+1,nfeval);
    for i = 1:nfeval
        f(:,i) = fun(x(:,i));
    end
end