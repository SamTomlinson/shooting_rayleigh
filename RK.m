%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 RK                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Uses a fourth order runge kutta method to march backward/forward in 
% time to get solution at opposite boudary



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
% bcs - values of boundary conditions (2D vector)
%
% a,b - two ends of the domain
% 
% k - spanwise wavenumber
% 
% eigval - eigenvalue shoot

% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                             Runge Kutta                               %



function [eta, p] = RK(a,b,deltaeta,bcs,rayleigh,baseU,...
        baseUdash,gamma,Tb,shoot1,c,beta)
    
    % Set up variables used in RK preallocate for speed
    
    eta = a:deltaeta:b; 
    n = length(eta);
    p = zeros(length(bcs),n);
    p(:,n) = bcs;
    
    % RK from far boundary in towards zero
    
    for k = n-1:-1:1
        j=5*k;
        k;
        
        k1 = rayleigh(eta(k+1),p(:,k+1),baseU(j-1),...
        baseUdash(j-1),gamma,Tb,shoot1,c,beta);
    
        k2 = rayleigh(eta(k+1)-0.5*deltaeta,p(:,k+1)-0.5*deltaeta*k1',...
            baseU(j-2),baseUdash(j-2),...
            gamma,Tb,shoot1,c,beta);
        
        k3 = rayleigh(eta(k+1)-0.5*deltaeta,p(:,k+1)-0.5*deltaeta*k2',...
            baseU(j-3),baseUdash(j-3),...
            gamma,Tb,shoot1,c,beta);
        
        k4 = rayleigh(eta(k+1)-deltaeta,p(:,k+1)-deltaeta*k3',...
            baseU(j-4),baseUdash(j-4),...
            gamma,Tb,shoot1,c,beta);
        
        p(:,k)=p(:,k+1)-deltaeta*(k1'+2*k2'+2*k3'+k4')/6;
    end
    
    
    