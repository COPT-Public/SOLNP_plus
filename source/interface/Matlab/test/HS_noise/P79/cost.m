function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
    %% P78
%     optx = [-1.717142,1.595708,1.827248,-.7636429,-.7636435]';
    fopt = .0787768209;
    f(1) = (x(1)-1)^2 + (x(1)-x(2))^2 + (x(2)-x(3))^2 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
    f(2) = x(1)+x(2)^2+x(3)^3-2-3*sqrt(2);
    f(3) = x(2) -x(3)^2 +x(4)+2-2*sqrt(2);
    f(4) = x(1)*x(5) -2;
    f = f';
    %
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
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && norm(f(2:length(f)))<tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
    f = f.*(1+1e-4*randn(4,1));
end