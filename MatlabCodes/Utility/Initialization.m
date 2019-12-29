

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization is a function for initializing the parameters of AGA1,
% AGA2, UGA1, UGA2, and NESUN. If some parameters specified by the user 
% Initialization uses these parameters. Otherwise, the default values 
% will be employed.
%
% INPUT: 
%
% x0                   % initial point
% options              % structure including the parameteres of schemes
%
%   .A                 % cell including matrices or linear operators
%   .MaxNumIter        % maximum number of iterations
%   .MaxNumFunEval     % maximum number of function evaluations
%   .MaxNumSubGradEval % maximum number of subgradient evaluations
%   .MaxNumLinOper     % maximum number of linear operators
%   .TimeLimit         % maximum running time
%   .epsilon           % accuracy parameter
%   .alpha_max         % maximum step size
%   .mu                % strong convexity parameter
%   .flag_x_error      % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
%   .flag_f_error      % 1 : saves f_error
%                      % 0 : do not saves f_error (default) 
%   .flag_time         % 1 : saves time for each iteration
%                      % 0 : do not saves time (default)  
%   .Stopping_Crit     % stopping criterion
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumSubGradEval is reached
%                      % 4 : stop if MaxNumLinOper is reached
%                      % 5 : stop if TimeLimit is reached
%                      % 6 : stop if eta <= epsilon
%                      % 7 : stop if Norm_dx/max(1,Norm_x) <= epsilon   
%
% OUTPUT:
%
% MaxNumIter           % maximum number of iterations
% MaxNumFunEval        % maximum number of function evaluations
% MaxNumSubGradEval    % maximum number of subgradient evaluations
% MaxNumLinOper        % maximum number of linear operators
% TimeLimit            % maximum running time
% epsilon              % accuracy parameter
% mu                   % strong convexity parameter
% x_opt                % optimizer 
% f_opt                % optimum 
% flag_x_error         % 1 : saves x_error
%                      % 0 : do not saves x_error (default)
% flag_f_error         % 1 : saves f_error
%                      % 0 : do not saves f_error (default)  
% flag_time            % 1 : saves time for each iteration
%                      % 0 : do not saves time (default)  
% Stopping_Crit        % stopping criterion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [epsilon,MaxNumIter,MaxNumFunEval,MaxNumLinOper, ...
    MaxNumSubGradEval,TimeLimit,flag_x_error,flag_f_error, ...
    flag_time,Stopping_Crit] = Initialization(options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Main body of Initialization.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(options,'epsilon') 
    epsilon = options.epsilon;
else
    epsilon = 10^(-5);
end

if isfield(options,'TimeLimit') 
    TimeLimit = options.TimeLimit;
else
    TimeLimit = inf;
end

if isfield(options,'MaxNumIter') 
    MaxNumIter = options.MaxNumIter;
else
    MaxNumIter = 5000;
end

if isfield(options,'MaxNumFunEval') 
    MaxNumFunEval = options.MaxNumFunEval;
else
    MaxNumFunEval = 10000;
end

if isfield(options,'MaxNumLinOper') 
    MaxNumLinOper = options.MaxNumLinOper;
else
    MaxNumLinOper = 15000;
end

if isfield(options,'MaxNumSubGradEval') 
    MaxNumSubGradEval = options.MaxNumSubGradEval;
else
    MaxNumSubGradEval = 10000;
end

if isfield(options,'flag_x_error') 
    flag_x_error = options.flag_x_error;
else
    flag_x_error = 0;
end

if isfield(options,'flag_f_error') 
    flag_f_error = options.flag_f_error;
else
    flag_f_error = 0;
end

if isfield(options,'flag_time') 
    flag_time = options.flag_time;
else
    flag_time = 0;
end

if isfield(options,'Stopping_Crit') 
    Stopping_Crit = options.Stopping_Crit;
else
    Stopping_Crit = 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of Initialization.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%