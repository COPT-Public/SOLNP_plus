% cnstr structure: 
%       pbl: lower bound for p
%       pbu: upper bound for p
%       ibl: lower bound of inequalities constrains
%       ibu: upper bound of inequalities constrains
%       p0:  initiate value of p
%       ib0: initiate value of inequality constrains
%
% op structure;
%       rho: penalty parameter
%       maxit: maximum outer iteration 
%       minit: maximum inner iteration
%       delta: norm of difference when caluculating gradient 
%       tol: stop condition 

function [p,jh,l,h,ic,count_cost]=solnp(cnstr,op,l,h)
%% setting initiate point and box constrains
if exist('cnstr','var')<=0.5
  error('Syntax error')
end
%
om=['SOLNP--> ';'         '];  %error message
% 
np = 0;
Ipc = [0,0];
Ipb=[0 0]; 
%lpb: Is the problem bounded
%np: dimension of x
if isfield(cnstr,'pbl')
    Ipc(1) = 1;
    Ipb(1) = 1;
    pb(:,1) = cnstr.pbl;
    [np,m] = size(cnstr.pbl);
    if m ~= 1
        error("Parameter array must one column.");
    end
end
if isfield(cnstr,'pbu')
    Ipc(2) = 1;
    Ipb(1) = 1;
    pb(:,2) = cnstr.pbu;
    [np,m] = size(cnstr.pbu);
    if m ~= 1
        error("Parameter array must one column.");
    end
end
if isfield(cnstr,'p0')
     p = cnstr.p0;
     [np,m] = size(cnstr.p0);
     if m ~= 1
        error("Parameter array must one column.");
     end
end
if Ipc(1) && Ipc(2)
    if ~isfield(cnstr,'p0')
        p = (pb(:,1) + pb(:,2))/2;
    end
    if max(p)==inf || min(p) == -inf
        error([om(1,:),'The user does not provide initiate point!']);
    end
else
    if ~isfield(cnstr,'p0')
        error([om(1,:),'The user does not provide initiate point!']);
    end
    if Ipc(1) == 0
        pb(:,1) = -inf * ones(np,1);
    end
    if Ipc(2) ==0
        pb(:,2) = inf * ones(np,1);
    end
end
%
if Ipb(1)>=0.5     %point is bounded
  if min(pb(:,2)-pb(:,1))<=0  %Error: lower bound is larger than upper bound
    error([om(1,:) 'The lower bounds of the parameter constraints ';...
          om(2,:) 'must be strictly less than the upper bounds.  ']);
  elseif min([p-pb(:,1);pb(:,2)-p])<=0 %Error: Infeasible initial point
    error([om(1,:) 'Initial parameter values must be within the bounds']);
  end
end
%% setting inequalities constrains
%ib: inequalities bound 
%ib0 estimated value of inequalities constrain at optimal point
%nic number of inquality constrains
ilc = 0; %is the problem lower constrained
iuc = 0; %is the problem upper constrained
nic = 0; %the number of inequality constrains
if isfield(cnstr,'ibl')
    ilc = 1;
    ib(:,1) = cnstr.ibl;
    [nic,m] = size(cnstr.ibl);
    if m ~= 1
        error("Inqualities array must one column.");
    end
end
if isfield(cnstr,'ibu')
    iuc = 1;
    ib(:,2) = cnstr.ibu;
    [nic,m] = size(cnstr.ibu);
    if m ~= 1
        error("Inqualities array must one column.");
    end
end
if isfield(cnstr,'ib0')
      ib0 = cnstr.ib0;
      [nic,m] = size(ib0);
      if m ~= 1
        error("Inqualities array must one column.");
      end
end
ic = ilc || iuc;

if ilc && iuc
    if ~isfield(cnstr,'ib0')
        ib0 = (ib(:,1)+ib(:,2))/2;
        if max(ib0) == inf || min(ib0) == -inf
            error([om(1,:),'The user does not provided initiate value of inequality constrains!']) ;
        end
    end
    if min(ib(:,2)-ib(:,1))<=0
      error([om(1,:) 'The lower bounds of the inequality constraints';...
            om(2,:) 'must be strictly less than the upper bounds. ']);
    end
elseif ic == 1
   if ~isfield(cnstr,'ib0')
        error([om(1,:),'The user does not provided initiate value of inequality constrains!']) ;
    end
   if ilc == 0
    ib(:,1) = -inf * ones(nic,1);
   end
   if iuc == 0 
    ib(:,2) = inf * ones(nic,1);
   end
end
if ic & min([ib0-ib(:,1);ib(:,2)-ib0])<=0
      error([om(1,:) 'Initial inequalities must be within the bounds']);
end


 %pb: inquality constrains and x constrains
 %p: estimated value of inquality constrains and x 
  if nic>=0.5
    if Ipb(1)>=0.5
      pb=[ib; pb];
    else
      pb=ib;
    end
    p=[ib0; p];
  end
clear cnstr
%lpb(2) is the problem bounded?(including point bound and function constrain)
if Ipb(1)+nic>=0.5
  Ipb(2)=1;
end
%% Setting paramters
% default optimization parameters
%replace default parameters if the user gives one
if ~exist('op','var'); op = struct();  end
if ~isfield(op,'rho'); op.rho = 1; end
if ~isfield(op,'maxit'); op.maxit = 10; end
if ~isfield(op,'minit'); op.minit = 10; end
if ~isfield(op,'delta'); op.delta = 1e-5; end
if ~isfield(op,'tol'); op.tol = 1e-4; end
if ~isfield(op,'tol_con'); op.tol_con = 1e-3; end
%if ~isfield(op,'random'); op.random = 0; end
rho=max(op.rho,0); maxit=op.maxit; minit=op.minit; delta=op.delta; tol=op.tol;
tolcon = op.tol_con;%random = op.random;
clear op
count_cost = 1;
%%  calculate initiate function value
ob=cost(p(nic+1:nic+np));
[m,n]=size(ob);
if n>1.5
  error([om(1,:) 'Cost function must return a column vector'])
end
if m<nic+1
  error([om(1,:) 'The number of constraints in your COST function does';...
        om(2,:) 'not match the number specified in the call to SOLNP.']);
end 
nec=m-1-nic;  %number of equality constrains
% if nec
%     Ipb(2) = 1 ;
% end
nc=m-1;      %number of total constrains
clear m n
%
j=ob(1);  
jh=j; t=0*ones(3,1);
%
if nc>0.5 
  if exist('l','var') <= 0.5 
    l=0*ones(nc,1); 
  end
  constraint=ob(2:nc+1);
  if nic>0.5
    if min([constraint(nec+1:nc)-pb(1:nic,1);...
            pb(1:nic,2)-constraint(nec+1:nc)]) > 0
      p(1:nic)=constraint(nec+1:nc);  %current point satisfy constrain
    end
    constraint(nec+1:nc)=constraint(nec+1:nc)-p(1:nic);
    % do not satisfy: difference with estimated value 
  end
  t(2)=norm(constraint);
  if max([t(2)-10*tol, nic])<=0
    rho=0;
  end 
else
  l=0;
end
%
if exist('h','var') <= 0.5 
  h=eye(np+nic);  %estiamted hessian
end
%% main iterations
mu=np; iteration=0;exit = 0;
while iteration<maxit
  iteration=iteration+1;

  op=[rho minit delta tol nec nic np Ipb tolcon];
  try
    [p,l,h,mu,count_cost]=subnp(p,op,l,ob,pb,h,mu,count_cost);
  catch ME
       disp(ME);
       exit = 1;
  end
%   disp('Updated parameters:')
%   p(nic+1:nic+np),
  ob=cost(p(nic+1:nic+np));
  count_cost = count_cost + 1;
  t(1)=(j-ob(1))/max(abs(ob(1)),1);
  j=ob(1);
  if nc>0.5
    constraint=ob(2:nc+1);
    if nic>0.5
      if min([constraint(nec+1:nc)-pb(1:nic,1);...
              pb(1:nic,2)-constraint(nec+1:nc)]) > 0
        p(1:nic)=constraint(nec+1:nc); 
      end
      constraint(nec+1:nc)=constraint(nec+1:nc)-p(1:nic);
    end
    t(3)=norm(constraint);
    if t(3)<10*tol
      rho=0; mu=min(mu, tol);
    end
    if t(3)<5*t(2) 
      rho=rho/5;
    elseif t(3)>10*t(2) 
      rho=5*max(rho, sqrt(tol)); 
    end
    if  max([tol+t(1), t(2)-t(3)]) <= 0 
      l=0*l; h=diag(diag(h)); 
    end 
    t(2)=t(3);
  end
  if exit == 1 || norm([t(1) t(2)]) <= tol 
    maxit=iteration; 
  end
  jh=[jh j];
end

%% output
if nic>0.5
  ic=p(1:nic);
end
p=p(nic+1:nic+np);
%
if norm([t(1) t(2)])<=tol
 fprintf([om(1,:) 'Completed in %d iterations\n'],iteration);
 fprintf([om(1,:),'The infeasibility is %e.\n'], t(3));
else
 fprintf([om(1,:) 'Exiting after maximum number of iterations. Tolerance not achieved\n']);
 fprintf([om(1,:),'The infeasibility is %e.\n'], t(3));
end
clear subnp_dual;
%
%-------------------------------------------------------------------------------
%  Variable Glossary:
%-------------------------------------------------------------------------------
%OB(1)      value of the cost objective function
%CONSTRAINT:vector of constraint values
%IB:        on input, contains the inequality constraint bounds + optionally
%           the values of the inequality constraints.  Gets converted to a 
%           NIC x 2 matrix of inequality constraint bounds.
%IC:        NIC x 1 vector of inequality constraints
%ITERATION: index for major iterations
%J:         previous value of the cost objective function
%JH:        history of the cost function
%LPB (2):   vector flag which indicates the presence of parameter bounds
%           and or inequality constraints.
%             LPB(1) refers to parameter bounds, it is 0 if there are 
%                    none, 1 if there are one or more.  
%             LPB(2) refers to constraints of either type.
%M:         number of rows of a matrix, usually temporary
%msg:       long error message string
%N:         number of columns of a matrix, usually temporary
%NC:        total number of constraints (=NEC+NIC)
%NEC:       number of equality constraints
%NIC:       number of inequality constraints
%NP:        number of parameters
%OPD:       vector of default optimization control variables
%om:        string for optimization messages
%OP:        vector of control variables.
%             It is passed in as:
%               [RHO MAJIT MINIT DELTA TOL] (all optional)
%             It is passed to ISISUBOPT as:
%               [RHO MAJIT MINIT DELTA TOL NEC NIC NP LPB(1) LPB(2)]
%P:         On input, contains the parameters to be optimized.  During the 
%           optimization this vector contains: [PIC;P]  where PIC contains 
%           pseudo parameters corresponding to the inequality constraints.
%PB:        on input, optionally contains the parameter bounds + optionally
%           the values of the parameters (one or the other or both can be 
%           specified).  Gets converted to a NPB x 2 matrix of parameter bounds.
%T (3):     vector of computed tolerances during optimization.
%             T(1) is the difference of the objective values between two 
%                  consecutive iterations
%             T(2) is NORM(CONSTRAINT) before a major iteration
%             T(3) is NORM(CONSTRAINT) after a major iteration
%
%-------------------------------------------------------------------------------
%  DOCUMENTATION:
%-------------------------------------------------------------------------------
%
%The Function SOLNP solves nonlinear programs in standard form:
%
%        minimize              J(P)
%        subject to            EC(P)  =0
%                   IB(:,1)<=  IC(P)  <=IB(:,2)
%                   PB(:,1)<=    P    <=PB(:,2).
%where
%
%  J       : Cost objective scalar function
%  EC      : Equality constraint vector function
%  IC      : Inequality constraint vector function
%  P       : Decision parameter vector
%  IB, PB  : lower and upper bounds for IC and P.
% 
% [P,L,JH,IC]=SOLNP(PB,IB,OP)
%
%  Output P      : optimal decision parameters
%         L      : optimal lagrangian multipliers
%         JH     : objective value's history of iterations
%         IC     : optimal values of inequality constraints
%        
%  Input  PB = [P {PB}]
%
%            P   : any initial parameter within the parameter bound
%            PB  : optional parameter bound
%                  PB(:,1) is the lower bound for P, and
%                  PB(:,2) is the upper bound for P.
%
%         IB = [{IC} IB] (optional)
%
%            IC  : (optional) best approximation values for inequality 
%                  constraints (IC) within the inequality constraint bound
%            IB  : inequality constraint bound
%                  IB(:,1) is the lower bound for IC, and
%                  IB(:,2) is the upper bound for IC.
%
%         OP=[RHO,MAJIT,MINIT,DELTA,TOL] (all optional with defaults)
%
%           RHO  : penalty parameter
%           MAJIT: maximum number of major iterations
%           MINIT: maximum number of minor iterations
%           DELTA: relative step size in forward difference evaluation
%           TOL  : tolerance on feasibility and optimality
%
%User-Defined Input function
%
%         OB : [OB]=COST(P,IT)
%
%           OB(1)        : cost objective function J value
%           OB(2:NEC+NIC): constraint function values of EC and IC
%           IT           : iteration number
%                           >0 - major iteration
%                           >0 - minor iteration
%                            0 - any other call