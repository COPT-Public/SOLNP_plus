clear cost
clear subnp_dual
clear subnp_reuse
clear subnp_sr1
clear
%cnstr.p0 = [-5,5,0]';
cnstr.pbl = [-4.5,-4.5,-5];
cnstr.pbu = -cnstr.pbl;
cnstr.ibl = 0;
cnstr.ib0 = 1;
op.tol = 1e-3;
op.tol_con = 1e-2;
op.min_iter = 3;
rep = 1;
info= SOLNP(cnstr,op);
cost(info.p,inf);