function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P106
    fopt = 680.6300573;
    f(1) = (x(1)-10)^2 + 5*(x(2)-12)^2 +x(3)^4 +3*(x(4)-11)^2 +10*x(5)^6 + 7*x(6)^2 +x(7)^4-4*x(6)*x(7)-10*x(6)-8*x(7);
    f(2) = 127 - 2*x(1)^2 - 3*x(2)^4 - x(3) -4*x(4)^2-5*x(5);
    f(3) = 282 - 7*x(1)-3*x(2)-10*x(3)^2 -x(4) +x(5);
    f(4) = 196 - 23*x(1) - x(2)^2 - 6*x(6)^2 + 8*x(7);
    f(5) = -4*x(1)^2 - x(2)^2 + 3*x(1)*x(2)-2*x(3)^2 - 5*x(6)+11*x(7);
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
   %f = f.*(1+1e-5*randn(7,1));
end