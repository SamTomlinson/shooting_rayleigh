%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           shooting_gotler                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Shooting method for solving 1D gotler stability equation with 
% Dirichlet BCS. 4th order RK and bisection to ensure BCs




%                                 Key                                 % 
%
% eta - grid points
%
% v - v0 and v0dash array 
%
% gotler - function containing de for gotler
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
% beta - streamwise wavenumber
%
% k - spanwise wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                         Output eigenfuntion                          %



function [eta, v,kval] = shooting_rayleigh2(rayleigh,deltaeta,tol,a,b,c)  

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; Tb=1; D=1; etab=1;
    B=-1.81447; 

     % Solve for the base flow 
    
    [~,baseT,baseTdash,baseU,baseUdash,baseUdashdash] ...
        = baseflow(C,Pr,D,etab,deltaeta,a,b);

    tic; % Begin time
    
    % Initialise vectors
    
    vec=[]; kvec=[];
    
    % Loop through different khat values 
    
    for shoot1=0.01:0.01:32

        % Far field boudary condition 
        
        a1 = [exp(-shoot1*b), -shoot1*exp(-shoot1*b)];
        
        % Runge kutta
        
        [~, F1] = RK(a,b,deltaeta,a1,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,c,gamma,Tb,shoot1);
    
        F1(1,1);
    
        % Boundary condition constraints
        
        H1=F1(2,1) + (-9*shoot1/a^4 + 2/a + ...
            B*shoot1/(a^(3-sqrt(7)))+ 2*C*shoot1)*F1(1,1);
        H2=F1(2,end) + shoot1*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[H1,vec];
        kvec=[shoot1,kvec];
    
    end
    
    % Plot H vs k
    
    figure('position', [0,0,800,800]); 
    plot(kvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{k})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Wavenumber, $\hat{k}$','Interpreter', 'LaTex','Fontsize',40)
    xlim([a,32])
    grid on
    hold off;
    
    % Calculate the corssing point
    vec;
    dif = abs(vec-0);
    match = dif == min(dif);
    idx = find(dif == min(dif));
    closest = vec(idx);
    kval= kvec(idx);
    
    % Plotting of solutions 
    
    a1 = [exp(-kval*b), -kval*exp(-kval*b)];
    [eta, F1] = RK(a,b,deltaeta,a1,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,c,gamma,Tb,kval);
     H1=F1(1,1)- exp(3*shoot1/a^3 + 2*log(a)...
     + B*shoot1*(-3-sqrt(7))/a^(4-sqrt(7))+ 2*C*shoot1);
     H2=F1(2,end) + shoot1*F1(1,end);
    
    %if (abs(H1) > tol)
    %    error('Zero bc not satiafied')
    %end
    %if (abs(H2) > tol)
    %    error('Far field bc not satisfied')
    %end
    
    v=F1;
    
    figure('position', [0,0,800,800]); 
    plot(eta,v(1,:),'k-','LineWidth',2);
    set(gca,'Fontsize',20)
    ylabel('Pressure. in the temp. adj. region $p_0$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([a,b])
    grid on
    hold off;
    
    toc
    
end