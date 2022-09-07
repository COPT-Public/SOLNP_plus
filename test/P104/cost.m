function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P106
    fopt = 3.9511634396;
    f(1) = 0.4*x(1)^(0.67)*x(7)^(-0.67)+0.4*x(2)^(.67)*x(8)^(-.67)+ 10 - x(1) - x(2);
    f(2) = 1 - 0.0588*x(5)*x(7) - 0.1*x(1);
    f(3) = 1 - 0.0588*x(6)*x(8) - 0.1*x(1)- 0.1*x(2);
    f(4) = 1 - 4*x(3)/x(5) - 2*x(3)^(-.71)/x(5)-0.0588*x(3)^(-1.3)*x(7);
    f(5) = 1 - 4*x(4)/x(6) - 2*x(4)^(-.71)/x(6)-0.0588*x(4)^(-1.3)*x(8);
    f(6) = 0.4*x(1)^(0.67)*x(7)^(-0.67)+0.4*x(2)^(.67)*x(8)^(-.67)+ 10 - x(1) - x(2);
   
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
   constraint = constraint + norm(min(x-0.1,0)+min(10-x,0))^2;
   constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
   %f = f.*(1+1e-5*randn(5,1));
end