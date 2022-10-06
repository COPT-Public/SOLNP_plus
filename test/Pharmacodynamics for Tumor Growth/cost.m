% This cost function calculate the dynamic system of tumor growth
% Input: [t_1,...,t_q,a_1,...,a_q]
% t_i: time of treatment
% a_i: amount of drug given at time t_i
% par indicator of where this function is called
% output: 
% sizeof tumor, maximum drug concentration, cumulative drug concentration.
function f = cost(x,par)
    persistent jh ch
    theta = [0.045, 4.52, 0.09, 0.11, 0.04, 0.00001, 0.09, 1];
    cmax = 0;
    cint = 0;

vmax = 1.1;
vcum = 65;
tmax = 200;
    res = [0,theta(7),theta(8),0];
    q = length(x)/2;
    %reorder the input 
    [t,idt] = sort(x(1:q));
    t = reshape(t,1,length(t));
    a = zeros(1,q);
    for i = 1:q
        a(i) = x(q+idt(i));
    end
    t = [0,t,tmax];
    a = [0,a];
    clear idt
    for i = 1:q+1
        %update C: give a(i) drug to the patient
        res(1) = res(1) + a(i);
        if t(i) == t(i+1)
            continue;
        end
        %calculate the C_max
        if cmax < res(1)
            cmax = res(1);
        end
        %calculate the integral if C
        cint = cint + res(1)*(1 - exp(-theta(1)*(t(i+1)-t(i))))/theta(1);
        %update P, Q, Q_P
        [tim,y] = ode45(@(time,treat)tumor(time,treat,t(i),res(1)),[t(i),t(i+1)],res(2:4)');
        %update initate point;
        res(2:4) = y(length(y(:,1)),:)';
        if exist("par","var")
            if par == inf
             plot(tim,sum(y,2));
            xlabel("Time");ylabel("Size of Tumor");
            hold on
%               subplot(1,2,1)
%               plotjh(jh,"Tumor Size");
%               hold on
%               subplot(1,2,2)
%               plotjh(ch,"Constraint violation");
%               hold on
          save("200-1.1-65 data SOLNP+.mat","jh","ch",'-mat');
            end
        end
        clear tim
        %update C(t)
        res(1) = res(1)*exp(-theta(1)*(t(i+1)-t(i)));
    end
    f = [sum(res(2:4));cmax;cint];
    if exist("par","var")
        if par == 1
            jh = [jh,f(1)];
        end
        if par == 2
            ch = [ch,sqrt((min(0,cmax)+min(vmax-cmax,0))^2+(min(0,cint)+min(vcum-cint,0))^2)];
        end
    else
        jh = [jh,f(1)];
        ch = [ch,sqrt((min(0,cmax)+min(vmax-cmax,0))^2+(min(0,cint)+min(vcum-cint,0))^2)];
    end
end

function f = tumor(t,x,ti,ci)
    P = x(1);
    Q = x(2);
    Qp = x(3);
    theta = [0.045, 4.52, 0.09, 0.11, 0.04, 0.00001, 0.09, 1];
    C = ci*exp(-theta(1)*(t - ti));
    K = 100;
    f(1) = theta(4)*P*(1 - (P+Q+Qp)/K) + theta(5)*Qp - theta(3)*P - theta(1)*theta(2)*C*P;
    f(2) = theta(3)*P-theta(1)*theta(2)*C*Q;
    f(3) = theta(1)*theta(2)*C*Q - theta(5)*Qp-theta(6)*Qp;
    f = f';
end

function plotjh(jh,yl)
    l = length(jh);
    plot(1:l,jh);
    xlabel("Function Evaluations");
    ylabel(yl);
    hold on
end