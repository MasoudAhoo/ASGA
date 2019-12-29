%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test_solvers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact 
clear
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nlist = {1};
Slist = {'ASGA-1','ASGA-2','ASGA-3','ASGA-4'};
flag  = 'ill_posed';
reg   = 'L1';
%reg   = 'L22L1';
%reg   = 'BCL1';
%reg   = 'BCL22L1';

Ttol          = 5;
maxit         = 1000;
Stopping_Crit = 5;
lambda        = 1e-2;
lambda1       = 1e-4;
lambda2       = 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1
    % ==================================================================
    if strcmp(flag,'ill_posed')
        % The function i_laplace.m is borrowed from the software package 
        % Regularization Tools
        n       = 1000; 
        [A,b,x] = i_laplace(n);
        xt      = x; 
        b       = b + 0.1*rand;
        %x0      = b;
        x0      = rand(n,1);
        
    end
    
    % =================== Start of implementations =====================
    switch reg 

      case 'L1'
        opt.b      = b;
        A          = {{A,A'},{}};
        opt.A      = A;
        opt.lambda = lambda;
          
        func1    = @ (varargin) L22L1R1(opt,varargin{:});
        subprob1 = @ (varargin) SubUnL1(lambda,varargin{:});
          
      case 'L22L1'
        opt.b       = b;
        A           = {{A,A'},{}};
        opt.A       = A;
        opt.lambda1 = lambda1;
        opt.lambda2 = lambda2;
        
        func1    = @ (varargin) L22L22L1R1(opt,varargin{:});
        subprob1 = @ (varargin) SubUnL1(lambda2,varargin{:});
        
      case 'BCL1'
        xl         = -ones(length(x0),1);
        xu         = ones(length(x0),1);
        opt.b      = b;
        A          = {{A,A'},{}};
        opt.A      = A;
        opt.lambda = lambda;
          
        func1    = @ (varargin) L22L1R1(opt,varargin{:});
        subprob1 = @ (varargin) SubBoL1(lambda,xl,xu,varargin{:});
         
      case 'BCL22L1'
        xl          = ones(length(x0),1);
        xu          = ones(length(x0),1);
        opt.b       = b;
        A           = {{A,A'},{}};
        opt.A       = A;
        opt.lambda1 = lambda1;
        opt.lambda2 = lambda2;
        
        func1    = @ (varargin) L22L22L1R1(opt,varargin{:});
        subprob1 = @ (varargin) SubBoL1(lambda2,xl,xu,varargin{:});
       
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Start of impelementations %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ======== determining initial guess for Lipschitz constant ========
    A1 = A{1}{1};
    if strcmp(reg,'L1')
        if strcmp(flag,'sparse')
            L = norm(A1,'fro')^2;
        else 
            L = norm(A1)^2;               % Lipschitz constant.
        end
    elseif strcmp(reg,'L22L1')
        if strcmp(flag,'sparse')
            L = norm(A1,'fro')^2+lambda1;
        else 
            L = norm(A1)^2+lambda1;       % Lipschitz constant.
        end
    elseif strcmp(reg,'BCL1')
        if strcmp(flag,'sparse')
            L = norm(A1,'fro')^2;
        else 
            L = norm(A1)^2;               % Lipschitz constant.
        end
    elseif strcmp(reg,'BCL22L1')
        if strcmp(flag,'sparse')
            L = norm(A1,'fro')^2+lambda1;
        else 
            L = norm(A1)^2+lambda1;       % Lipschitz constant.
        end
    end
        
    [fx0 gx0] = func1(x0);
    z         = 0.5*x0;              
    [fz gz]   = func1(z);
    L0        = norm(gx0-gz)/norm(x0-z);
    
    options.L                 = L; 
    options.L0                = L0; 
    options.A                 = A;
    options.MaxNumIter        = maxit;
    options.MaxNumFunEval     = 1000;
    options.MaxNumLinOper     = 300;
    options.MaxNumSubGradEval = 5000;
    options.Stopping_Crit     = Stopping_Crit;
    options.TimeLimit         = Ttol;
    options.epsilon           = 1e-1; 
    if strcmp(reg,'L1') 
        options.mu = 0;
    elseif strcmp(reg,'L22L1')
        options.mu = lambda1;
    elseif strcmp(reg,'BCL1') 
        options.mu = 0;
        options.xl = xl;
        options.xu = xu;
    elseif strcmp(reg,'BCL22L1')
        options.mu = lambda1;
        options.xl = xl;
        options.xu = xu;
    end
                
    % ======================== applying solvers ========================
    for j = 1:length(Slist) 

        switch Slist{j}               
                
            case 'ASGA-1'
                options.nu  = 1;
                options.Lnu = L;
                Lk_finder=@(varargin)Holder_constant_finder(varargin{:});
                fprintf('Running ASGA-1 ...\n')
                tic
                [x,f,out] = ASGA1( func1,Lk_finder,subprob1,x0,options );         
                toc
                F = out.F;
                F_best = F(1);
                for k = 2 : length(F)
                    if F(k) < F_best
                        F_best = F(k);
                    end
                    F(k)   = F_best;
                end
                fun_asga1   = F(end);
                t_asga1     = out.T;
                nfunc_asga1 = out.Nfunc;
                f_eval{j}   = F';
            
            case 'ASGA-2'
                fprintf('Running AAGA-2 ...\n')
                tic
                [ x,f,out ] = ASGA2( func1,subprob1,x0,options );         
                toc
                F = out.F;
                F_best = F(1);
                for k = 2 : length(F)
                    if F(k) < F_best
                        F_best = F(k);
                    end
                    F(k)   = F_best;
                end
                fun_asga2   = F(end);
                t_asga2     = out.T;
                nfunc_asga2 = out.Nfunc;
                f_eval{j}   = F';
            
            case 'ASGA-3'
                options.nu  = 1;
                options.Lnu = L;
                Lk_finder=@(varargin)Holder_constant_finder(varargin{:});
                fprintf('Running ASGA-3 ...\n')
                tic
                [ x,f,out ]=ASGA3(func1,Lk_finder,subprob1,x0,options);         
                toc
                F = out.F;
                F_best = F(1);
                for k = 2 : length(F)
                    if F(k) < F_best
                        F_best = F(k);
                    end
                    F(k)   = F_best;
                end
                fun_asga3   = F(end);
                t_asga3     = out.T;
                nfunc_asga3 = out.Nfunc;
                f_eval{j}   = F';
            
            case 'ASGA-4'
                fprintf('Running ASGA-4 ...\n')
                tic
                [ x,f,out ] = ASGA4( func1,subprob1,x0,options );         
                toc
                F = out.F;
                F_best = F(1);
                for k = 2 : length(F)
                    if F(k) < F_best
                        F_best = F(k);
                    end
                    F(k)   = F_best;
                end
                fun_asga4   = F(end);
                t_asga4     = out.T;
                nfunc_asga4 = out.Nfunc;
                f_eval{j}   = F';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun_asga1
fun_asga2
fun_asga3
fun_asga4

nfunc_asga1
nfunc_asga2
nfunc_asga3
nfunc_asga4


figure(1)

for k = 1: length(Slist)
    switch Slist{k}
                  
        case 'ASGA-1'
            xk = floor(exp((1:20)*(log(length(f_eval{k}))/20)));
            Fk = f_eval{k}(xk);
            loglog(xk,Fk,'vk','LineWidth',2)

        case 'ASGA-2'
            xk = floor(exp((1:20)*(log(length(f_eval{k}))/20)));
            Fk = f_eval{k}(xk);
            loglog(xk,Fk,'oc','LineWidth',2)
            
        case 'ASGA-3'
            xk = floor(exp((1:15)*(log(length(f_eval{k}))/15)));
            Fk = f_eval{k}(xk);
            loglog(xk,Fk,'*r','LineWidth',2)

        case 'ASGA-4'
            xk = floor(exp((1:15)*(log(length(f_eval{k}))/15)));
            Fk = f_eval{k}(xk);
            loglog(xk,Fk,'pb','LineWidth',2)
    end
    hold on
    xlabel('iterations');
    ylabel('function values');
    legend('ASGA-1','ASGA-2','ASGA-3','ASGA-4');
end
for k = 1: length(Slist)
    switch Slist{k}
                  
        case 'ASGA-1'
            loglog(f_eval{k},'-k','LineWidth',1.7)

        case 'ASGA-2'
            loglog(f_eval{k},'-c','LineWidth',1.7)
            
        case 'ASGA-3'
            loglog(f_eval{k},'-r','LineWidth',1.7)

        case 'ASGA-4'
            loglog(f_eval{k},'-b','LineWidth',1.7)
    end
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of the script %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

