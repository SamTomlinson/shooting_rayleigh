%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               rayleigh                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Encode gotler system




%                                 Key                                 % 
%
% eta - grid points
%
% q - derivatives p_0 and p_0'
%
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE
%
% khat - streamwise wavenumber
%
% c - eigen value 
%
% gamma, Tb - flow parameters
%



%                          Rayleigh system                             %

function vecq = rayleigh(~,q,baseU,...
        baseUdash,gamma,Tb,shoot1,c,beta)

% Diff of q1 is q2 

vecq(1) = q(2);
    
% Diff of q2 is the rest of the system
    
vecq(2) = 2*(0.5*(gamma-1)*(Tb+1)*baseUdash)*q(2)/...
        ((0.5*(gamma-1)*(Tb+1)*baseU)-c) ...
        + ((shoot1^2+beta^2)^2)*((1-(0.5*(gamma-1)*(Tb+1)*baseU))^2)*q(1);
    
    
    
    
    