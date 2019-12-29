

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SubUnL1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SubUnL1: this function solves the following problem 
%
%          min    f(y)+<g_f(y),x-y>+L/2||x-y||_2^2+lambda ||x||_1, 
%          s.t.   x in V,
%
%          which appears in Nesterov-type optimal gradient methods for 
%          convex composite objective functions.
%
% INPUT:
%
% y                    % point y in the above minimization problem
% lambda               % regularization parameter
% mu                   % strong convexity parameter
% optsub               % structure including the parameteres
%
%   .Sum_si            % Scaling parameter for the subproblem 2
%   .Sum_sigi          % sum of scaling parameters and gradients
%                      % for the subproblem 2
%
% OUTPUT:
%
% x                    % the best approximation of the minimizer          
%            
% WRITTEN BY: 
%
% Masoud Ahookhosh
% Faculty of Mathematics, University of Vienna, Austria
%
% LAST UPDATE: 
%
% April 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = SubUnL1( lambda,mu,y,optsub )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Error messages for input and output %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5
    error('The number of input arguments is more than what is needed');
elseif nargin < 4
    error('The number of input arguments is not enough');
end;

if isempty(y)
    error('SubUnL1 needs the y to be defined');
elseif ~isa(y,'numeric')
    error('y should be a numeric vector');
end

if isempty(lambda)
    error('SubUnL1 needs the lambda to be defined');
elseif ~isa(lambda,'numeric')
    error('lambda should be a scalar');
end

if isempty(mu)
    error('SubUnL1 needs the mu to be defined');
elseif ~isa(mu,'numeric')
    error('mu should be a scalar');
end

if isempty(flag)
    error('SubUnL1 needs flag (1 or 2) to be defined');
elseif ~isa(flag,'numeric')
    error('flag should be a scalar');
end

if isempty(optsub)
    error('SubUnL1 needs the structure array optsub to be defined');
elseif ~isa(optsub,'struct')
    error('optsub should be a structure array');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Main body of SubUnL1.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if flag == 1
%     % This section solves the primal composite gradient subproblem
%     
%     sk1      = optsub.sk1;
%     grad   = optsub.g;
%     
%     tau = sk1;
%     n   = length(y);
%     x   = zeros(n,1);
%     for i = 1 : n
%         if y(i) > tau*(grad(i)+lambda)
%             x(i) = y(i)-tau*(grad(i)+lambda);
%         elseif y(i) < tau*(grad(i)-lambda)
%             x(i) = y(i) - tau*(grad(i)-lambda);
%         end
%     end
%        
% elseif flag == 2
%     % This section solves the dual composite gradient subproblem
%     
%     Sum_si = optsub.Sum_si;
%     tau    = lambda * Sum_si;
%     grad   = optsub.Sum_sigi;
%     
%     n   = length(y);
%     x   = zeros(n,1);
%     for i = 1 : n
%         if y(i) > grad(i) + tau
%             x(i) = y(i) - grad(i) - tau;
%         elseif y(i) < grad(i) - tau
%             x(i) = y(i) - grad(i) + tau;
%         end
%     end
%     
% else
%     error('flag is not correct')
% end

Sum_si   = optsub.Sum_si;
Sum_sigi = optsub.Sum_sigi;

% n   = length(y);
% x   = zeros(n,1);
% c1  = y-Sum_sigi-lambda*Sum_si;
% c2  = y-Sum_sigi+lambda*Sum_si;
% 
% ind1 = find(c1>0);
% ind2 = find(c2<0);
% 
% x(ind2) = c2(ind2);
% x(ind1) = c1(ind1);

%y = prox_{tau*|.|_1}(x) = max(0,1-tau/|x|)*x
%
%   Proximal operator for the scalar L1 norm.
%
%   Copyright (c) 2010 Gabriel Peyre
r   = 1+mu*Sum_si;
yt  = (y-Sum_sigi)/r;
tau = (lambda*Sum_si)/r;

x = max(0,1-tau./max(abs(yt),1e-10)).*yt;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of SubUnL1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%