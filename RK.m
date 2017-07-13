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



function [eta, p] = RK(a,b,deltaeta,bcs,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,gamma,Tb,shoot1,khat)
    
    % Set up variables used in RK preallocate for speed
    
    eta = a:deltaeta:b; 
    n = length(eta);
    p = zeros(length(bcs),n);
    p(:,n) = bcs;
    
    % RK from far boundary in towards zero
    
    for k = n-1:-1:1
        
        k;
        
        k1 = rayleigh(eta(k+1),p(:,k+1),baseT(k+1),baseTdash(k+1),baseU(k+1),...
        baseUdash(k+1),baseUdashdash(k+1),gamma,Tb,shoot1,khat);
    
        k2 = rayleigh(eta(k+1)-0.5*deltaeta,p(:,k+1)-0.5*deltaeta*k1',baseT(k+1),...
            baseTdash(k+1),baseU(k+1),baseUdash(k+1),baseUdashdash(k+1),...
            gamma,Tb,shoot1,khat);
        
        k3 = rayleigh(eta(k+1)-0.5*deltaeta,p(:,k+1)-0.5*deltaeta*k2',baseT(k+1),...
            baseTdash(k+1),baseU(k+1),baseUdash(k+1),baseUdashdash(k+1),...
            gamma,Tb,shoot1,khat);
        
        k4 = rayleigh(eta(k+1)-deltaeta,p(:,k+1)-deltaeta*k3',baseT(k+1),baseTdash(k+1),...
            baseU(k+1),baseUdash(k+1),baseUdashdash(k+1),...
            gamma,Tb,shoot1,khat);
        p(:,k)=p(:,k+1)-deltaeta*(k1'+2*k2'+2*k3'+k4')/6;
    end
    
    
    