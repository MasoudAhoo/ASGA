

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% StopCriterion.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StopCriterion is a function checking that one of the stopping criteria 
% holds to terminate AGA1, AGA2, UGA1, UGA2, and NESUN. It also perepare 
% the status determining why the algorithm is stopped.
%
% INPUT:
% 
% fx                    % function value at the current point x
% Niter                 % number of iterations 
% Nfunc                 % number of function evaluations
% Nsubgrad              % number of subgradient evaluations
% Nlinop                % number of linear operators
% Norm_dx               % norm(xb-x_old)
% Time                  % running time
% MaxNumIter            % maximum number of iterations
% MaxNumFunEval         % maximum number of function evaluations
% MaxNumSubGradEval     % maximum number of subgradient evaluations
% MaxNumLinOper         % maximum number of linear operators
% TimeLimit             % maximum running time
% epsilon               % accuracy parameter
% Stopping_Crit         % stopping criterion
%
%                       % 1 : stop if MaxNumIter is reached
%                       % 2 : stop if MaxNumFunEval is reached
%                       % 3 : stop if MaxNumSubGradEval is reached
%                       % 4 : stop if MaxNumLinOper is reached
%                       % 5 : stop if TimeLimit is reached
%                       % 6 : stop if (fx-f_min) <= epsilon
%                       % 7 : stop if Norm_dx/max(1,Norm_x) <= epsilon   
%
% OUTPUT:
%
% StopFlag              % 1: if one of the stopping criteria holds
%                       % 0: if none of the stopping criteria holds
% Status                % the reason of termination
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function [StopFlag, Status] = StopCriterion(fx,Niter,Nfunc,Nlinop, ...
    Nsubgrad,Norm_dx,Norm_x,Time,MaxNumIter,MaxNumFunEval, ...
    MaxNumLinOper,MaxNumSubGradEval,TimeLimit,epsilon,Stopping_Crit)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Main body of StopCriterion.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Stopping_Crit
  case 1
    if Niter >= MaxNumIter
      StopFlag = 1;
      Status   = 'Maximum number of iterations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 2
    if Nfunc >= MaxNumFunEval
      StopFlag = 1;
      Status   = 'Maximum number of function evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 3
    if Nlinop >= MaxNumLinOper
      StopFlag = 1;
      Status   = 'Maximum number of linear operators is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 4
    if Nsubgrad >= MaxNumSubGradEval
      StopFlag = 1;
      Status   = 'Maximum number of subgradient evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 5
    if Time >= TimeLimit
      StopFlag = 1;
      Status   = 'Time limit is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 6
    if (fx-f_min) <= epsilon
      StopFlag = 1;
      Status   = '|f(x_k) - f_min| <= epsilon';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 7
    if flag && Norm_dx/max(1,Norm_x) <= epsilon
      StopFlag = 1;
      Status = '||x_k - x_k-1|| / max(1,norm(x_k)) <= epsilon';
    else
      StopFlag = 0;
      Status   = [];
    end

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of StopCriterion.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
