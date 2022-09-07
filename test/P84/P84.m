clear cost
clear
prob.p0 =  [2.52,2,37.5,9.25,6.8];% + 1e-3*randn(5,1);
prob.ibu = [294000,294000,277200]';
prob.ibl = [0,0,0]';
prob.pbu = [1000,2.4,60,9.3,7]';
prob.pbl = [0,1.2,20,9,6.5]';
op.max_iter = 50;
op.min_iter = 10;
op.tol = 1e-3;
op.tol_con = 1e-2;
op.tol_restart = 1;
op.re_time = 5;
info= SOLNP(prob,op);
cost(info.p,inf);