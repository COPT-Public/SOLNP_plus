% This program use NOMAD to solve the tumor growth problem
% Input: initiate time and drug dosage
% Ouput: optimal time and drug dosage
clear cost
nq = 4;
vmax = 1.1;
vcum = 65;
tmax = 200;
x0 = ones(1,nq)*tmax/2;
x0 = [x0,0.5*ones(1,nq)]';
lb = zeros(2*nq,1);
ub = [tmax*ones(nq,1);ones(nq,1)];
param = struct('display_degree','0','bb_output_type','OBJ PB PB PB PB','max_bb_eval','3000','min_poll_size','1e-8');
f = @(x)fun(x,vmax,vcum);
tic
[x0,fval,hinf,exit_status,nfeval] = nomad(f,x0,lb,ub,param);
toc
cost(x0,inf)

function f = fun(x,vmax,vcum)
    f = cost(x);
    f = [f(1);-f(2,1);-f(3,1);f(2:3,1)-[vmax;vcum]];
%     f = [f(1); max(f(2:3))];
%     f = f'; 
end
