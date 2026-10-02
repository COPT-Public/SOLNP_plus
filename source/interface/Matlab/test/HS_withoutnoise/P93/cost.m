function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
%     persistent jh;
    tol = 1e-2;
  %% P93
    fopt = 135.075961;
    f(1) = .0204*x(1)*x(4)*(x(1)+x(2)+x(3)) + .0187*x(2)*x(3)*(x(1)+1.57*x(2)+x(4))...
        +.0607*x(1)*x(4)*x(5)^2*(x(1)+x(2)+x(3))+.0437*x(2)*x(3)*x(6)^2*(x(1)+1.57*x(2)+x(4));
    f(2) = .001*x(1)*x(2)*x(3)*x(4)*x(5)*x(6) - 2.07;
    f(3) = 1 - .00062*x(1)*x(4)*x(5)^2*(x(1)+x(2)+x(3))-.00058*x(2)*x(3)*x(6)^2*(x(1)+1.57*x(2)+x(4));
    f = f';

    if exist('par','var') && par==inf
        fprintf('COST-->  Final result: F(x) = %e\n',f(1));
        fprintf('COST-->  The algorithm has called function for %d times\n',calltimes);
%         plot(1:calltimes,jh);
%         hold on
        if finaltime ~= inf
            fprintf('COST-->  The algorithm requires %d evaluations to get an optimal solution\n\n ',finaltime);
        else
            fprintf('COST-->  Fail to get an optimal solution\n');
        end
        return;
    end
%     jh = [jh;f(1)];
    if isempty(calltimes)
           calltimes = 0;
           tag = 1;
    end
    calltimes = calltimes+1;
    if isempty(finaltime)
           finaltime = inf;
    end
   constraint = norm(min(f(2:length(f)),zeros(2,1)))^2;
   constraint = constraint + norm(min(x,0));
   constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
%      f = f.*(1+1e-4*randn(3,1));
end