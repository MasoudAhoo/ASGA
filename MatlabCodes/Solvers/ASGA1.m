

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASGA1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ x,f,out ] = ASGA1(func,Lk_finder,subprob,x0,options) 
% ASGA1 is an accelerated gradient algorithm for solving convex  
%      composite optimization
%                          min f(x) + psi(x)  s.t.  x in C,
%      where both f and phi are convex functions and C is a convex set.
%
% INPUT:
%
% func                 % function handle for the objective function
% Lk_finder            % a routine finding Lk
% subprob              % function handle for associated subproblems
% x0                   % initial point
% options              % structure including the parameteres of scheme
%
%   .A                 % cell including matrices or linear operators
%   .MaxNumIter        % maximum number of iterations
%   .MaxNumFunEval     % maximum number of function evaluations
%   .MaxNumSubGradEval % maximum number of subgradient evaluations
%   .MaxNumLinOper     % maximum number of linear operators
%   .TimeLimit         % maximum running time
%   .epsilon           % accuracy parameter
%   .mu                % strong convexity parameter
%   .L                 % Holder constant
%   .nu                % level of smoothness
%   .x_opt             % optimizer 
%   .f_opt             % optimum
%   .flag_x_error      % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
%   .flag_f_error      % 1 : saves f_error
%                      % 0 : do not saves f_error (default)
%   .flag_time         % 1 : saves f_error
%                      % 0 : do not saves f_error (default)
%   .Stopping_Crit     % stopping criterion
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumSubGradEval is reached
%                      % 4 : stop if MaxNumLinOper is reached
%                      % 5 : stop if TimeLimit is reached
%                      % 6 : stop if (fx-f_min) <= epsilon
%                      % 7 : stop if Norm_dx/max(1,Norm_x) <= epsilon
%
% OUTPUT:
%
% x                    % the best approximation of the optimizer
% f                    % the best approximation of the optimum
% out                  % structure including more output information
%
%   .T                 % running time
%   .Niter             % total number of iterations
%   .Nfunc             % total number of function evaluations
%   .Nsubgrad          % total number of subgradient evaluations
%   .Nlinop            % total number of employed linear operators
%   .F                 % array including all function values            
%   .x_error           % relative error norm(xb(:)-x_opt(:))/norm(x_opt)
%   .f_error           % relative error (fk-fs)/(f0-fs))    
%   .Status            % reason of termination
%
% REFERENCE: 
%
% [1] M. Ahookhosh, Accelerated first-order methods for large-scale 
%     convex minimization, Submitted, (2016)
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


function [ x,f,out ] = ASGA1( func,Lk_finder,subprob,x0,options )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Initializing and setting the parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

% ================ Error messages for input and output =================
if nargin > 5
    error('The number of input arguments is more than what is needed');
elseif nargin < 5
    error('The number of input arguments is not enough');
end;

if isempty(func)
    error('ASGA1 needs the function handle func to be defined');
elseif ~isa(func,'function_handle')
    error('func should be a function handle');
end

if isempty(subprob)
    error('ASGA1 needs the function handle subprob to be defined');
elseif ~isa(subprob,'function_handle')
    error('subprob should be a function handle');
end

if isempty(x0)
    error('ASGA1 needs an starting point x0 to be defined');
elseif ~isa(x0,'numeric')
    error('x0 should be a numeric vector');
end

% =================== initializing the parameters ======================
% ===== user has requested viewing the default values of "options" =====
[epsilon,MaxNumIter,MaxNumFunEval,MaxNumLinOper,MaxNumSubGradEval, ...
    TimeLimit,flag_x_error,flag_f_error,flag_time,Stopping_Crit] ...
    = Initialization(options);

if isfield(options,'x_opt')
    x_opt=options.x_opt;
elseif flag_x_error==1
    error('x_error requires to x_opt be specified');
end

if isfield(options,'nu')
    nu=options.nu;
else 
    error('nu should be be specified');
end

if isfield(options,'Lnu')
    Lnu = options.Lnu;
else 
    error('Lnu should be be specified');
end

if isfield(options,'mu')
    mu = options.mu;
else 
    error('mu should be be specified');
end

if flag_x_error == 1
    Nxopt      = sqrt(sum(x_opt(:).^2));
    x_error(1) = sqrt(sum((x0(:)-x_opt(:)).^2))/Nxopt;
end

if flag_f_error == 1
    f_error(1) = 1;
end

if flag_time == 1
    Time(1) = 0;
end

xk       = x0;
zk       = x0;
Niter    = 1;
[fx0]    = func(x0);
Nfunc    = 1;
Nlinop   = 1;    
Nsubgrad = 1;
F        = fx0;
Sum_skgk = 0;
StopFlag = 0;
Sk       = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of ASGA1.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;

% ======================= start of the main loop =======================
while ~StopFlag
      
    % ================= finding Lk and computing sk1 ===================
    r   = 1+Sk*mu;
    if nu == 0 
        Lkt = (1/(2*r*epsilon))*Lnu^2;
        Lk = 2*(1+2*Sk)*Lkt;
    elseif nu == 1
        Lk = Lnu;
    else
        if Niter == 1
            L0 = Lnu;
        else
            L0 = Lk;
        end
        Lk  = Lk_finder(Sk,r,nu,Lnu,epsilon,L0);  
    end
    sk1 = (r+sqrt(r^2+4*Lk*Sk*r))/(2*Lk);
       
    % ====== computing yk and a (sub)gradient of f ======
    Sk1       = Sk+sk1;
    alphak    = sk1/Sk1;
    yk        = alphak*zk+(1-alphak)*xk;
    [fyk gyk] = func(yk);
    Nfunc     = Nfunc+1;
    Nlinop    = Nfunc+2;     
    Nsubgrad  = Nfunc+1;
    
    % ====================== computing zk1 and xk1 ====-================
    Sum_skgk        = Sum_skgk+sk1*(gyk-mu*yk); 
    optsub.Sum_si   = Sk;
    optsub.Sum_sigi = Sum_skgk;
    zk1             = subprob(mu,x0,optsub);  
    xk1             = alphak*zk1+(1-alphak)*xk;
    fxk1            = func(xk1);
    
    % ====================== updating parameters =======================
    xk_old   = xk;
    xk       = xk1;
    Niter    = Niter+1
    fxk      = fxk1;
    F(Niter) = fxk;    
    zk       = zk1;
    Sk       = Sk1;
      
    % ================== checking stopping criteria ====================
    Time     = toc(T0);
    if Stopping_Crit == 7
        Norm_dxk = norm(xk-xk_old);
        Norm_xk  = norm(xk);
    else
        Norm_dxk = 0;
        Norm_xk  = 0;
    end
                
    [StopFlag, Status] = ...
         StopCriterion(fxk,Niter,Nfunc,Nlinop,Nsubgrad,Norm_dxk, ...
         Norm_xk,Time,MaxNumIter,MaxNumFunEval,MaxNumLinOper, ...
         MaxNumSubGradEval,TimeLimit,epsilon,Stopping_Crit);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status
x            = xk;
f            = fxk;
T            = Time;
out.T        = T;
out.F        = F';
out.Niter    = Niter;
out.Nfunc    = Nfunc;
out.Nsubgrad = Nsubgrad;
out.Nlinop   = Nlinop;
out.Status   = Status;
        
if flag_x_error == 1
    out.x_error = x_error;          
end
if flag_f_error == 1
    out.f_error = f_error;          
end
if flag_time == 1
    out.Time = Time;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of ASGA1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%