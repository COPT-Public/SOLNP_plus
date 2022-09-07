function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;

%% P28
optx = [log(log(10)),log(10),10];
fopt = -log(log(10));
f = [-x(1),x(2)-exp(x(1)),x(3)-exp(x(2))];
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
    constraint = norm(min(f(2:length(f)),zeros(length(f )-1,1)))^2;
    constraint = constraint + norm(min(x,[0,0,0]')+min([100,100,10]'-x,0))^2;
    constraint = sqrt(constraint);
   if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
end