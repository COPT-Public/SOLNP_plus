clear cost
clear
rng(223)
prob.p0 = [-2,1]';
op.tol = 1e-4;
op.tol_con = 1e-4;
op.min_iter = 1;
op.maxit = 50;

fun.cost = @(x,nfeval)cost_para(x,nfeval,1);
fun.grad = @grad;
info= SOLNP(prob,op,fun);

function f = cost_para(x,nfeval,nc)
    f = zeros(nc+1,nfeval);
    for i = 1:nfeval
        f(:,i) = cost(x(:,i));
    end
end
function gradf = gradf(x)
    gradf = [2*(x(1)-1);0];
end
function gradg = gradg(x)
    gradg = 10*[-2*x(1);1];
end
function g = grad(x)
    g = [gradf(x);gradg(x)];
end