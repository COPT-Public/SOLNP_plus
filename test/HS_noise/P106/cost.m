function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P106
    fopt = 7049.330923;
    f(1) = x(1)+x(2)+x(3);
    f(2) = 1 - 0.0025*(x(4)+x(6));
    f(3) = 1 - 0.0025*(x(5)+x(7) - x(4));
    f(4) = 1 - 0.01*(x(8)-x(5));
    f(5) = x(1)*x(6) - 833.33252*x(4) -100*x(1) +83333.333;
    f(6) = x(2)*x(7) - 1250*x(5) -x(2)*x(4) +1250*x(4);
    f(7) = x(3)*x(8) - 1250000 - x(3)*x(5) +2500*x(5);
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
   constraint = norm(min(f(2:length(f)),0))^2;
%    constraint = constraint + norm(min(x,0));
   constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
   f = f.*(1+1e-4*randn(7,1));
end