function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
   %% P83
         optx = [78,33,29.99526,45,36.77581]';
         fopt = -30665.53867;
         a = [85.334407,.0056858,.0006262, .0022053, 80.51249 ,.0071317,.0029955,.0021813,9.300961,.0047026,.0012547,.0019085];
         f = [5.3578547*x(3)^2 + .8356891*x(1) * x(5) + 37.293239*x(1)- 40792.141 , a(1) + a(2)*x(2)*x(5) + a(3)*x(1)*x(4) - a(4)*x(3)*x(5),a(5) + a(6)*x(2)*x(5) + a(7)*x(1)*x(2) + a(8)*x(3)^2 - 90,a(9) + a(10)*x(3)*x(5) + a(11)*x(1)*x(3) + a(12)*x(3)*x(4) - 20];

 f = f';
    if exist('par','var') && par==inf
        fprintf('COST-->  Final result: F(x) = %e\n',f(1));
        fprintf('COST-->  The algorithm has called function for %d times\n',calltimes);
        if finaltime ~= inf
            fprintf('COST-->  The algorithm requires %d evaluations to get an optimal solution\n\n ',finaltime);
        else
            fprintf('COST-->  Fail to get an optimal solution\n');
        end
        return;
    end
    if isempty(calltimes)
           calltimes = 0;
           tag = 1;
    end
    calltimes = calltimes+1;
    if isempty(finaltime)
           finaltime = inf;
    end
    constraint = norm(min(f(2:length(f)),[0;0;0])+ max(f(2:length(f))-[92;20;5],zeros(3,1)))^2;
    constraint = constraint + norm(min(x-[78,33,27,27,27]',0)+min([102,45,45,45,45]'-x,0));
    constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
%     f = f.*(1 + 1e-4*randn(4,1));
   % f(1) = f(1)*(1+1e-4*randn(1,1));
end