function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P38
    fopt = 0;
    f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1 - x(3))^2 + 10.1*((x(2)-1)^2 + (x(4)-1)^2) + 19.8*(x(2)-1)*(x(4)-1);
    f = f';
    if isempty(finaltime)
        
        finaltime = inf;
    end
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
 
    constraint = norm(min(x + 10 * ones(4,1),0')+min(10 * ones(4,1) -x,0));
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
     noise = 1e-4*randn(1,1);
%     disp(noise);
    f = f.*(1+noise);
end