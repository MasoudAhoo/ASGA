

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Holder_constant_finder.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Holder_constant_finder is a function computing the solution the 
% nonlineat equation
%
% L-(1+Sk*mu+((1+Sk*mu)^2+4*LkbSk*(1+Sk*mu))^(1/2))^((1-nu)/(1+nu))Lt=0
%
% INPUT:
% 
% Sk      % the parameter Sk
% mu      % strong convexity parameter
% nu      % level of smoothness
% Lnu     % Holder constant
% epsilon % accuracy parameter
% alphak  % step-size
% L0      % initial guess for zero finder (the last Lk)
%
% OUTPUT:
%
% Lk       % the parapeter Lk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Lk = Holder_constant_finder(Sk,r,nu,Lnu,epsilon,L0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Main body of Holder_constant_finder.m %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 6
    error('The number of input arguments is inconsistent');
end

Lkt = ((1-nu)/(2*r*epsilon*(1+nu)))^((1-nu)/(1+nu))*Lnu^(2/(1+nu));

fun = @(L) L-(r+sqrt(r^2+4*L*Sk*r))^((1-nu)/(1+nu))*Lkt;
Lk = fzero(fun,L0);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End of Holder_constant_finder.m %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

