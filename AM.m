%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 MS                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Uses a fourth order multistep method to march backward/forward in 
% time to get solution at opposite boudary



%                                 Key                                 % 
%
% eta - grid points
%
% p - p0 and v0dash array 
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



function [eta, p] = AM(a,b,deltaeta,bcs1,bcs2,bcs3,bcs4,rayleigh,baseU,...
        baseUdash,gamma,Tb,shoot1,c,beta)
    
    % Set up variables used in RK preallocate for speed
    
    eta = [a:deltaeta:b,b+deltaeta,b+2*deltaeta,b+3*deltaeta,b+4*deltaeta]; 
    n = length(eta);
    p = zeros(length(bcs1),n);
    p(:,n-3) = bcs1;
    p(:,n-2) = bcs2;
    p(:,n-1) = bcs3;
    p(:,n) = bcs4;
    
    % Additional points for base flow
    
    baseU=[baseU,0,0,0,0]; baseUdash=[baseUdash,0,0,0,0];
    
    
    % RK from far boundary in towards zero
    
    for k = n-4:-1:1
        k;
        k+1;
        k1 = rayleigh(eta(k+1),p(:,k+1),baseU(k+1),...
        baseUdash(k+1),gamma,Tb,shoot1,c,beta);
        k+2;
        k2 = rayleigh(eta(k+2),p(:,k+2),...
            baseU(k+2),baseUdash(k+2),...
            gamma,Tb,shoot1,c,beta);
        k+3;
        k3 = rayleigh(eta(k+3),p(:,k+3),...
            baseU(k+3),baseUdash(k+3),...
            gamma,Tb,shoot1,c,beta);
        k+4;
        k4 = rayleigh(eta(k+4),p(:,k+4),...
            baseU(k+4),baseUdash(k+4),...
            gamma,Tb,shoot1,c,beta);
        
        p(:,k)=p(:,k+1)-deltaeta*((55/24)*k1'-(59/24)*k2'...
            +(37/24)*k3'-(3/8)*k4');
    end