function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P60
    fopt = .9535288567;
    f(1) = (x(1)-x(2))^2 +(x(1)+x(2)-10)^2/9 + (x(3)-5)^2;
    f(2) = 48 - x(1)^2 - x(2)^2 -x(3)^2;
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

    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && norm(min(f(2),0))<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
 %f = f.*(1+1e-4*randn(2,1));
end