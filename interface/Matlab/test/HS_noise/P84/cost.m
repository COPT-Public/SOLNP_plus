function f = cost(x,par)    
    persistent calltimes;
    persistent finaltime;
    persistent tag;
    tol = 1e-2;
  %% P84
    optx = [4.53743097 , 2.4 , 60 , 9.3 , 7];
    fopt = -5280335.133;
    a = [-24345,-8720288.849,150512.5253,-156.6950325,476470.3222,729482.8271,-145421.402,2931.1506,-40.427932,5106.192,15711.36,-155011.1084,4360.53352,12.9492344,10236.884,13176.786,-326669.5104,7390.68412,-27.8986976,16643.076,30988.146];
    f = -a(1)-a(2)*x(1)-a(3)*x(1)*x(2)-a(4)*x(1)*x(3)-a(5)*x(1)*x(4)-a(6)*x(1)*x(5);
    f(2) = a(7)*x(1)+a(8)*x(1)*x(2)+a(9)*x(1)*x(3)+a(10)*x(1)*x(4)+a(11)*x(1)*x(5);
    f(3) = a(12)*x(1)+a(13)*x(1)*x(2)+a(14)*x(1)*x(3)+a(15)*x(1)*x(4)+a(16)*x(1)*x(5);
    f(4) = a(17)*x(1)+a(18)*x(1)*x(2)+a(19)*x(1)*x(3)+a(20)*x(1)*x(4)+a(21) *x(1)*x(5);

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
   constraint = norm(min(f(2:length(f)),[0;0;0])+ max(f(2:length(f))-[294000;294000;277200],zeros(3,1)))^2;
    constraint = constraint + norm(min(x -[0,1.2,20,9,6.5]' ,0)+min([1000,2.4,60,9.3,7]'-x,0))^2;
    constraint = sqrt(constraint);
    if abs(f(1)-fopt)/max(abs(fopt),1) <= tol && constraint<=tol  && tag == 1
        finaltime = calltimes;
        tag = 0;
    end
      f = f.*(1+ 1e-4*randn(4,1));
end