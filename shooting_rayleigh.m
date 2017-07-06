%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           shooting_rayleigh                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Shooting method for solving 1D rayleigh stability equation with 
% Dirichlet BCS. 4th order RK and bisection to ensure BCs



%                                 Key                                 % 
%
% eta - grid points
%
% p - p and pdash array 
%
% rayleigh - function containing de for gotler
%
% deltaeta - step size
%
% tol - tolerance
%
% bcs - values of boundary conditions (2D vector)
%
% init - initial guess, if not specified given as [-5,10]
%
% a,b - two ends of the domain
%
% flow parameters - gamma (specific heat), Pr (prandtl), C (sutherlands
% constant), D (fitting parameter), etab (matching point for edge of
% adjustment region)
%
% khat - combined wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                               Example                                %
%
% [eta, v] = shooting_rayleigh(@rayleigh,h,zero,a,b,con,type,init,khat)
% [eta, v] = shooting_rayleigh(@rayleigh,0.006,1e-6,1,7,[0,0],'df',1);
%
% i.e. solve bvp in gotler on [1,7] with bcs v(1)=0 and v(7)=0, 
% tolerance 1e-6, wavenumber =1 khat=0.1




%                         Output eigenfuntion                          %



function [eta, p] = shooting_rayleigh(rayleigh,deltaeta,tol,a,b,bcs,...
    init,khat) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; Tb=1; D=1; etab=1; M=1; c=-0.993937;
    kappa=0.1;

     % Solve for the base flow 
    
    [~,baseT,baseTdash,baseU,baseUdash,baseUdashdash] ...
        = baseflow(C,Pr,D,etab,deltaeta,a,b);

    tic; % Begin time
    
    % If my number of arguements is 8 then initial guesses gave been 
    % specified if not take these to be -1 and 1.
    
    if nargin == 19
        shoot1 = init(1); shoot2 = init(2);
    else
        shoot1 = -5; shoot2 = 10;
    end
    
    % Sets up boundary condition vectors
    
    a1 = [shoot1 bcs(1)];
    a2 = [shoot2 bcs(1)]; 
    
   % Now iterate solution outwards using Rk method 
    
    [~, F1] = RK(a,b,deltaeta,a1,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,c,gamma,Tb,khat); 
    [~, F2] = RK(a,b,deltaeta,a2,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,c,gamma,Tb,khat);       
    
     F1 = F1(1,end) - bcs(2);
     F2 = F2(1,end) - bcs(2);   
    
    % Identify if as root is possible by checking for sign change
    
    % Check 
    % F1*F2
    
    if (F1*F2 > 0) 
        error('The root does not exist')
    end
    
    % Set one shoot for iteration 
    
    F3 = F1;
   
    % Iteration to home in on axis crossing
    
    while (abs(F3) > tol) 
        
        % Check
        % F3
        
        % Bring one shoot in half the distance between the teo
        
        shoot3 = (shoot1 + shoot2)/2;

        % Renforce conditions and rerun RK solver on loop adjusting 
        % to compensate for average overshooting root
        
        a3 = [bcs(1) shoot3];                      
        
        [eta, F3] = RK(a,b,deltaeta,a3,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,c,gamma,Tb,khat); 
        
        % Check
        % F3(r,end);
        
        p = F3; F3 = F3(1,end) - bcs(2); 
        if (F1*F3 < 0)
            shoot2 = shoot3; F2 = F3;            
        elseif (F1*F3 > 0)
            shoot1 = shoot3; F1 = F3;
        else
            error('Something has gone horribly wrong, probs NANS');           
        end
        
    end      
    
    % Plotting of solutions 
    
    figure('position', [0,0,800,800]); 
    plot(eta,p(1,:),'k-','LineWidth',2); hold on; 
    plot(eta,p(2,:),'r-','LineWidth',2); 
    set(gca,'Fontsize',20)
    l1=legend('$p(\eta)$','$p_{\eta}(\eta)$');
    set(l1, 'Interpreter','LaTex','Fontsize',30);
    ylabel('Pres. in the temp. adj. region $p$','Interpreter', 'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([1,7])
    grid on
    hold off;
    toc
end