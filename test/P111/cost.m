function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P111
    fopt = -47.76109026;
    c = [-6.089,-17.164,-34.054,-5.914,-24.721,-14.986,-24.100,-10.708,-26.662,-22.179];
    f(1) = 0;
    for i=1:10
        f(1) = f(1) +  exp(x(i))*(c(i) + x(i) - log(sum(exp(x))));
    end
    f(2) = exp(x(1)) + 2*exp(x(2))+ 2*exp(x(3))+ exp(x(6))+ exp(x(10))-2;
    f(3) = exp(x(4)) + 2*exp(x(5)) +exp(x(6))+exp(x(7))-1;
    f(4) = exp(x(3)) + exp(x(7)) + exp(x(8))+ 2*exp(x(9))+ exp(x(10))-1;
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
   constraint = norm(f(2:length(f)))^2;
%    constraint = constraint + norm(min(x,0));
   constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
   %f = f.*(1+1e-5*randn(7,1));
end