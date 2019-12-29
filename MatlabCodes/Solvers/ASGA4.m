

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASGA4.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ x,f,out ] = ASGA4(func,subprob,x0,options) 
% ASGA4 is an accelerated gradient algorithm for solving convex  
%      composite minimization
%                      min f(x) + psi(x)  s.t.  x in C,
%      where both f and phi are convex functions and C is a convex set.
%
% INPUT:
%
% func                 % function handle for the objective function
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
%   .L0                % initial approximation of Holder constant
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
%                      % 7 : stop if Norm_dx/max(1,Norm_x)<=epsilon
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


function [ x,f,out ] = ASGA4( func,subprob,x0,options )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Initializing and setting the parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

% ================ error messages for input and output =================
if nargin > 4
    error('The number of input arguments is more than what is needed');
elseif nargin < 3
    error('The number of input arguments is not enough');
end;

if isempty(func)
    error('ASGA4 needs the function handle func to be defined');
elseif ~isa(func,'function_handle')
    error('func should be a function handle');
end

if isempty(subprob)
    error('ASGA4 needs the function handle subprob to be defined');
elseif ~isa(subprob,'function_handle')
    error('subprob should be a function handle');
end

if isempty(x0)
    error('ASGA4 needs an starting point x0 to be defined');
elseif ~isa(x0,'numeric')
    error('x0 should be a numeric vector');
end

% =================== initializing the parameters ======================
% ===== user has requested viewing the default values of "options" =====
[epsilon,MaxNumIter,MaxNumFunEval,MaxNumLinOper,MaxNumSubGradEval, ...
    TimeLimit,flag_x_error,flag_f_error,flag_time,Stopping_Crit] ...
    = Initialization(options);

if isfield(options,'x_opt')
    x_opt = options.x_opt;
elseif flag_x_error == 1
    error('x_error requires to x_opt be specified');
end

if isfield(options,'L0')
    L0 = options.L0;
else 
    error('L0 should be be specified');
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

Lk       = L0;
yk       = x0;
vk       = x0;
Niter    = 1;
fy0      = func(yk);
Nfunc    = 1;
Nlinop   = 1;    
Nsubgrad = 1;
F        = fy0;
Sk       = 0;
Sum_skgk = 0;
StopFlag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of ASGA4.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0 = tic;

% ============================ outer loop ==============================
while ~StopFlag
      
    % ========================= computing sk1 ==========================
    Lkb = Lk;
    r   = 1+Sk*mu;   
    sk1 = (r+sqrt(r^2+4*Lkb*Sk*r))/(2*Lkb);
       
    % ============= computing yk and a (sub)gradient of f ==============
    Sk1         = Sk+sk1;
    alphak      = sk1/Sk1;
    xk1         = alphak*vk+(1-alphak)*yk;
    [fxk1 gxk1] = func(xk1);
    Nfunc       = Nfunc+1;
    Nlinop      = Nlinop+2;    
    Nsubgrad    = Nsubgrad+1;
    
    % ====================== computing uk1 and yk1 ====-================
    optsub.Sum_si   = sk1;
    optsub.Sum_sigi = sk1*(gxk1-mu*xk1);
    uk1             = subprob(mu,vk,optsub);  
    yk1             = alphak*uk1+(1-alphak)*yk;
    fyk1            = func(yk1);
    
    % ========================== inner loop ============================
    k=0;
    yk1_xk1 = yk1-xk1;
    while fyk1>fxk1+gxk1'*yk1_xk1+0.5*Lkb*sum(yk1_xk1.^2)+ ...
          0.5*epsilon*alphak

        if k > 10
           break;
        else
            k=k+1;
        end

        % ======================= computing sk1 ========================
        Lkb = 4*Lkb;
        %Lkb = 3*Lkb;
        r   = 1+Sk*mu;   
        sk1 = (r+sqrt(r^2+4*Lkb*Sk*r))/(2*Lkb);
           
        % ==== computing the sequence y and gradient of smooth part ====
        Sk1         = Sk+sk1;
        alphak      = sk1/Sk1;
        xk1         = alphak*vk+(1-alphak)*yk;
        [fxk1 gxk1] = func(xk1);
        Nfunc       = Nfunc+1;
        Nlinop      = Nlinop+2;    
        Nsubgrad    = Nsubgrad+1;
    
        % ======= solving the auxiliary problem and computing xk =======
        optsub.Sum_si   = sk1;
        optsub.Sum_sigi = sk1*(gxk1-mu*xk1);
        uk1             = subprob(mu,vk,optsub);    
        yk1             = alphak*uk1+(1-alphak)*yk;   
        fyk1            = func(yk1);
        yk1_xk1         = yk1-xk1;
        
    end
    % ======================= end of inner loop ========================
    
    % ====================== updating parameters =======================
    yk_old   = yk;
    yk       = yk1;
    Niter    = Niter+1
    fyk      = fyk1;
    F(Niter) = fyk;      
    Sk       = Sk1;
    Lk       = 0.6*Lkb; 
    
    % ========================= computing vk1 ==========================  
    Sum_skgk        = Sum_skgk+sk1*(gxk1-mu*xk1); 
    optsub.Sum_si   = Sk;
    optsub.Sum_sigi = Sum_skgk;
    vk              = subprob(mu,x0,optsub);
       
    % ================== checking stopping criteria ====================
    Time = toc(T0);
    if Stopping_Crit == 7
        Norm_dyk = norm(yk-yk_old);
        Norm_yk  = norm(yk);
    else
        Norm_dyk = 0;
        Norm_yk  = 0;
    end
                
    [StopFlag, Status] = ...
         StopCriterion(fyk,Niter,Nfunc,Nlinop,Nsubgrad,Norm_dyk, ...
         Norm_yk,Time,MaxNumIter,MaxNumFunEval,MaxNumLinOper, ...
         MaxNumSubGradEval,TimeLimit,epsilon,Stopping_Crit);
    
end
% ========================= end of outer loop ==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status
x            = yk;
f            = fyk;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of ASGA4.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%