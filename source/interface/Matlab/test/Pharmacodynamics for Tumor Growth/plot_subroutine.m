
solnp_lsty = "b--";
pdfo_lsty = "r-";
nomad_lsty = "g:";

% subplot(1,2,1)
load("200-1.1-65 data PDFO.mat");
plotjh(jh,ch,"Tumor Size",1,pdfo_lsty);
hold on
load("200-1.1-65 data SOLNP+.mat");
plotjh(jh,ch,"Tumor Size",1,solnp_lsty);
load("200-1.1-65 data NOMAD.mat");
plotjh(jh,ch,"Tumor Size",1,nomad_lsty);
legend("PDFO","SOLNP+","NOMAD");


% subplot(1,2,2)
% load("200-1.1-65 data PDFO.mat");
% plotjh(ch,"Infeasibility",1,pdfo);
% hold on
% load("200-1.1-65 data SOLNP.mat");
% plotjh(ch,"Infeasibility",1,solnp);
% load("200-1.1-65 data NOMAD.mat");
% plotjh(ch,"Infeasibility",1,nomad);
% legend("PDFO","DFNP","NOMAD");


function plotjh(jh,ch,yl,start,co)
   
    l = length(jh);
    tol = 1e-1;
    jh_best = zeros(l,1);
    %select the best solution within the 
    if ch(1) > tol
        jh_best(1) = 20;
    else
        jh_best(1) = jh(1);
    end
    for i=2:l
        if ch(i) > tol
            jh_best(i) = jh_best(i-1);
        else
            jh_best(i) = min(jh_best(i-1),jh(i));
        end
    end
    jh_best = [jh_best;jh_best(end)];

    plot([start:l,3000],jh_best(start:l+1),co,"LineWidth",1.2);
    xlabel("Function Evaluations");
    ylabel(yl);
    hold on
end