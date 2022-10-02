function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;

%% P46
optx = ones(5,1);
fopt = 0;
f = [(x(1)-x(2))^2+ (x(3)-1)^2 + (x(4)-1)^4+ (x(5)-1)^6 ,x(1)*x(4)*x(1)+sin(x(4)-x(5))-1,x(2)+x(4)^2*x(3)^4-2];

 f = f';
 %
 %f = f + 1e-4*randn(3,1);
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
    if abs(f-fopt)/max(abs(fopt),1) <= tol & norm(f(2:length(f)))<tol & tag == 1
        finaltime = calltimes;
        tag = 0;
       % disp('Reach!');
    end
%         f = f.*(1 + 1e-4*randn(3,1));
end