
solnp = "b--";
pdfo = "r-";
nomad = "g:";

% subplot(1,2,1)
% load("200-1.1-65 data PDFO.mat");
% plotjh(jh,"Tumor Size",1000,pdfo);
% hold on
% load("200-1.1-65 data SOLNP+.mat");
% plotjh(jh,"Tumor Size",1000,solnp);
% load("200-1.1-65 data NOMAD.mat");
% plotjh(jh,"Tumor Size",1000,nomad);
% legend("SOLNP+","NOMAD");


% subplot(1,2,2)
load("200-1.1-65 data PDFO.mat");
plotjh(ch,"Infeasibility",1,pdfo);
hold on
load("200-1.1-65 data SOLNP+.mat");
plotjh(ch,"Infeasibility",1,solnp);
load("200-1.1-65 data NOMAD.mat");
plotjh(ch,"Infeasibility",1,nomad);
legend("PDFO","SOLNP+","NOMAD");


function plotjh(jh,yl,start,co)
   
    l = length(jh);
    plot(start:l,jh(start:l),co,"LineWidth",1.2);
    xlabel("Function Evaluations");
    ylabel(yl);
    hold on
end