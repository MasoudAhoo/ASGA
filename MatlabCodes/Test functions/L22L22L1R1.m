

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L22L22L1R1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L22L22L1R1 is a function generating both function values and   
% (sub)gradient evaluations of the convex test function:
% 
%           f(x) = 1/2||Ax-b||_2^2+1/2||x||_2^2+lambda*||x||_1.
%
% INPUT:
%
% xk         % current point;
% opt        % structure includes required parameters;
%    .b      % observation;
%
% OUTPUT:
%
% fk % function value of f at xk
% gk % (sub)gradient of f at xk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fk gk] = L22L22L1R1(opt,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main body of L22L22L1R1.m %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 2
    error('The number of input arguments is not valid');
end   

if nargout >= 3 
    error('The number of output arguments is not valid');
end

xk  = varargin{1};
   
% ======================== Function evaluation =========================
b       = opt.b;
A       = opt.A;
lambda1 = opt.lambda1;
lambda2 = opt.lambda2;

Axk   = feval(LOPMVM(A),xk,1);
Axk_b = Axk{1}{1}-b;
fk    = 0.5*sum(Axk_b.^2)+0.5*lambda1*sum(xk.^2)+lambda2*norm(xk,1);
    
% ===================== (Sub)gradient evaluation =======================
if nargout > 1
    gk1 = {{Axk_b},{}};
    gk2 = {lambda1*xk,{}};
    gk  = SubGradEval(A,gk1,gk2);
end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of L22L22L1R1.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%