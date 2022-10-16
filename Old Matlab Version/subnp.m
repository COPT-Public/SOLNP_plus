function [p,y,h,l,count_cost]=subnp(p0,op,yy,ob,pb,h,l,count_cost)
%
rho=op(1); maxit=op(2); delta=op(3); tol=op(4);   nec=op(5); 
nic=op(6); np=op(7);    nc=nec+nic;  npic=np+nic; lpb=op(8:9); ch=1;
clear op
alp=[0 0 0];
%
% make the scale for the cost, the equality constraints, the inequality
% constraints, and the parameters
%
if nec>0   
  scale=[ob(1);ones(nec,1)*max(abs(ob(2:nec+1)))];
else        
  scale=1;
end
if lpb(2)<=0 
  scale=[scale; p0];
else        
  scale=[scale; ones(size(p0))];
end
% scale=[scale; ones(size(p0))];
scale=min(max(abs(scale),tol),1/tol);

% scale the cost, the equality constraints, the inequality constraints, 
% the parameters (inequality parameters AND actual parameters), 
% and the parameter bounds if there are any
% Also make sure the parameters are no larger than (1-tol) times their bounds
%
ob=ob./scale(1:nc+1);  
p0=p0./scale(nec+2:nc+np+1);
if lpb(2)>=0.5
  if lpb(1)<=0.5
    mm=nic;
  else 
    mm=npic;
  end
  pb=pb./[scale(nec+2:nec+mm+1) scale(nec+2:nec+mm+1)];
end
%
% scale the lagrange multipliers and the Hessian
%
if nc>0.5
  yy=scale(2:nc+1).*yy/scale(1);
end
h=h.*(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)')/scale(1);
%
om=['SOLNP--> ';'         '];
msg=[om(1,:) 'Redundant constraints were found. Poor              '
     om(2,:) 'intermediate results may result.  Suggest that you  '
     om(2,:) 'remove redundant constraints and re-OPTIMIZE.       '];
%
j=ob(1);
a=[0*ones(nec,nic);-eye(nic)]; % Transform the inequatlity into eqality constrains
g=0*ones(npic,1);
%  g: estimated gradient using difference
%  a: Jacobi Matrix of ALL Constrain functions
if nc>0.5
  constraint=ob(2:nc+1);
  for i=1:np
    p0(nic+i)=p0(nic+i)+delta;
    ob=cost(p0(nic+1:npic).*scale(nc+2:nc+np+1))./scale(1:nc+1);
    count_cost = count_cost + 1;
    g(nic+i)=(ob(1)-j)/delta;
    a(:,nic+i)=(ob(2:nc+1)-constraint)/delta;
    p0(nic+i)=p0(nic+i)-delta;
  end
  if nic>0.5
    constraint(nec+1:nec+nic)=constraint(nec+1:nec+nic)-p0(1:nic);
    % constrain: value of all equality constrains
  end
  if cond(a)>1/eps
    error(msg);
  end
  b=a*p0-constraint;
end
%% find a feasible solution
msg=[om(1,:) 'The linearized problem has no feasible     '
     om(2,:) 'solution.  The problem may not be feasible.'];
%
if nc>0.5 % if there are any constrains
  ch=-1; % ch use to record whether the concurrent point is feasible
  alp(1)=tol-max(abs(constraint)); %how much constrains are not obeyed 
  if alp(1)<=0
    ch=1;
    if lpb(2)==0
      p0=p0-a'*((a*a')\constraint);
      alp(1)=1;
    end
  end
  if alp(1)<=0 %current point is different from estimated value
    p0(npic+1)=1;  
    a=[a, -constraint];
    c=[zeros(1,npic), 1];
    dx=ones(npic+1,1); 
    go=1; %flag whether 
    minit=0;
    while go>=tol
      % Find a feasible solution of linearized constrains 
      % mm: total number of inequalities constrains including equations and
      % box constrains
      minit=minit+1;
      gap=[p0(1:mm,1)-pb(:,1),pb(:,2)-p0(1:mm,1)];
      gap=sort(gap,2);
      dx(1:mm)= min(gap(:,1),10^4); %lower bound 
      dx(npic+1)=p0(npic+1); % dx(npic+1)?? what 
      if lpb(1)==0 %no box constrains
        dx(mm+1:npic)=max([dx(1:mm);100])*ones(npic-mm,1);
      end
      y=(a*diag(dx))'\(dx.*c'); 
      v=dx.*(dx.*(c'-a'*y)); %fuction of v??
      if v(npic+1)>0
        z=p0(npic+1)/v(npic+1);
        for i=1:mm
          if v(i)<0
            z=min(z,-(pb(i,2)-p0(i))/v(i));
          elseif v(i)>0
            z=min(z,(p0(i)-pb(i,1))/v(i)); 
          end
        end
        if z>= p0(npic+1)/v(npic+1)
          p0=p0-z*v;   % 
        else
          p0=p0-0.9*z*v; 
        end
        go=p0(npic+1);
        if minit >= 10 %reach maximum iteration
          go=0; 
        end
      else
        go=0;
        minit=10;
      end
    end
    if minit>=10
      disp(msg),
    end
    a=a(:,1:npic); 
    b=a*p0(1:npic);
  end
end
%
clear constraint c z v gap;
%
p=p0(1:npic);
y=0; 
if ch>0
  ob=cost(p(nic+1:npic).*scale(nc+2:nc+np+1))./scale(1:nc+1);
  count_cost = count_cost +1;
end
j=ob(1);
%
if nic>0.5
  ob(nec+2:nc+1)=ob(nec+2:nc+1)-p(1:nic);
end
if nc>0.5
  ob(2:nc+1)=ob(2:nc+1)-a*p+b; %why? b = a*p??
  j=ob(1)-yy'*ob(2:nc+1)+rho*norm(ob(2:nc+1))^2; % augmented lagrangian yy: multiplier
end
%
minit=0;
%% main iterations
while minit<maxit %Use SQP to solve subproblem
  minit=minit+1;
  if ch>0 % different with estimation
    for i=1:np %caluculate gradient (again) using difference of function value at p
               %(We obtain a new fe asible point now )
      p(nic+i)=p(nic+i)+delta;
      obm=cost(p(nic+1:npic).*scale(nc+2:nc+np+1))./scale(1:nc+1);
      count_cost = count_cost +1;
      if nic>0
        obm(nec+2:nc+1)=obm(nec+2:nc+1)-p(1:nic);
      end
      if nc>0
        obm(2:nc+1)=obm(2:nc+1)-a*p+b;
        obm=obm(1)-yy'*obm(2:nc+1)+rho*norm(obm(2:nc+1))^2;
      end
      g(nic+i)=(obm-j)/delta; %g: gradient
      p(nic+i)=p(nic+i)-delta;
    end
    if nic>0.5
      g(1:nic)=0*yy(nec+1:nc); %gradient variables related to the ineqalities are zero
    end
  end
  if minit>1 %start from the second iteration: Update the approximated Hessian
    yg=g-yg; %yg difference of gradients
    sx=p-sx; %difference of x
    sc(1)=sx'*h*sx;
    sc(2)=sx'*yg;
    if sc(1)*sc(2)>0
      sx=h*sx;
      h=h-(sx*sx')/sc(1)+(yg*yg')/sc(2); % BFGS update of hessian
    end
  end
  dx=0.01*ones(npic,1); 
  if lpb(2)>0.5
    gap=[p(1:mm,1)-pb(:,1),pb(:,2)-p(1:mm,1)];
    gap=sort(gap,2);% Why sort?
    gap=gap(:,1)+(1e-8)*ones(mm,1);
    dx(1:mm,1)=max(ones(mm,1)./gap,0.01);
    if lpb(1)<=0
      dx(mm+1:npic)=min([dx(1:mm);0.01])*ones(npic-mm,1);
    end
  end
  go=-1;
  l=l/10;

  while go<=0 %use interior ellipsold algothrim to solve QP
    c=chol(h+l*diag(dx.*dx)); % h+l*diag(dx.*dx) may not be spd?
    c=inv(c);
    yg=c'*g;
    if nc<=0
      u=-c*yg;
    else
      y=(c'*a')\yg;
      u=-c*(yg-(c'*a')*y);
    end
    p0=u(1:npic)+p;
    if lpb(2)<=0
      go=1;
    else
      go=min([p0(1:mm)-pb(:,1);pb(:,2)-p0(1:mm)]);
      l=3*l; % adjust the radius
    end
  end
  %y = zeros(nec,1);
  %% line search
  alp(1)=0;ob1=ob;ob2=ob1;sob(1)=j;sob(2)=j;
  pt(:,1:2)=[p p];alp(3)=1.0;pt(:,3)=p0;
  ob3=cost(pt(nic+1:npic,3).*scale(nc+2:nc+np+1))./scale(1:nc+1);
  count_cost = count_cost +1;
  sob(3)=ob3(1);
  if nic>0.5
    ob3(nec+2:nc+1)=ob3(nec+2:nc+1)-pt(1:nic,3);
  end
  if nc>0.5
    ob3(2:nc+1)=ob3(2:nc+1)-a*pt(:,3)+b;
    sob(3)=ob3(1)-yy'*ob3(2:nc+1)+rho*norm(ob3(2:nc+1))^2;
  end
  %unorm = norm(u);
  go=1;
  lstime = 0 ;
  
  while go>tol 
    lstime = lstime +1;
    alp(2)=(alp(1)+alp(3))/2;
    pt(:,2)=(1-alp(2))*p+alp(2)*p0;
    ob2=cost(pt(nic+1:npic,2).*scale(nc+2:nc+np+1))./scale(1:nc+1);
    count_cost = count_cost +1;
    sob(2)=ob2(1);
    if nic>0.5
      ob2(nec+2:nc+1)=ob2(nec+2:nc+1)-pt(1:nic,2);
    end
    if nc>0.5
      ob2(2:nc+1)=ob2(2:nc+1)-a*pt(:,2)+b;
      sob(2)=ob2(1)-yy'*ob2(2:nc+1)+rho*norm(ob2(2:nc+1))^2;
    end
    obm=max(sob);
    if obm<j
      obn=min(sob);
      go=tol*(obm-obn)/(j-obm);
    end
    if sob(2)>=sob(1)
      sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
    elseif sob(1)<=sob(3)
      sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
    else
      sob(1)=sob(2);ob1=ob2;alp(1)=alp(2);pt(:,1)=pt(:,2);
    end
    if go>=tol
      go=alp(3)-alp(1);
    end
  end
  sx=p;yg=g;ch=1;
  obn=min(sob); % new minimum value
  if j<=obn  %j: old minimum value
    maxit=minit;
  end
  reduce=(j-obn)/(1+abs(j));  
  if reduce<tol
    maxit=minit;
  end
  if sob(1)<sob(2)
    j=sob(1);p=pt(:,1);ob=ob1;
  elseif sob(3)<sob(2)
    j=sob(3);p=pt(:,3);ob=ob3;
  else 
    j=sob(2);p=pt(:,2);ob=ob2;
  end
  clear ob1 ob2 ob3 pt;
end
p=p.*scale(nec+2:nc+np+1);  % unscale the parameter vector
if nc>0.5
  y=scale(1)*y./scale(2:nc+1); % unscale the lagrange multipliers
end
h=scale(1)*h./(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)');
%
if reduce>tol
  disp([...
  om(1,:) 'Minor optimization routine did not converge in the ';...
  om(2,:) 'specified number of minor iterations.  You may need';...
  om(2,:) 'to increase the number of minor iterations.        '])
end
end

%% gradient function

%
%------------------------------------------------------------------------------
%  Variable Glossary:
%------------------------------------------------------------------------------
%
%ALP (4):   vector of ??. Alp(4) is "alpha"
%CH:        "change"
%CONSTRAINT:vector of constraint values
%G:         gradient
%GAP (2):   lower and upper "gap"
%LPB (2):   vector flag which indicates the presence of parameter bounds
%           and or inequality constraints.
%             LPB(1) refers to parameter bounds, it is 0 if there are 
%                    none, 1 if there are one or more.  
%             LPB(2) refers to constraints of either type.
%msg:       long error message string
%NC:        total number of constraints (=NEC+NIC)
%NEC:       number of equality constraints
%NIC:       number of inequality constraints
%NP:        number of parameters
%om:        string for optimization messages
%P0:        parameter vector on input
%P:         parameter vector on output
%SOB (3):   vector of ??
%Z (5):     vector of ??