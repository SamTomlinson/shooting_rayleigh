%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 loop                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% loops improving accuracy of shoot



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
% adjustment region), Tb (wall temp), c (neutral phase speed)
%
% khat - combined wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE
%

function [eigval H1] = loop(eigval,a,b,deltaeta,rayleigh,baseT,...
    baseTdash,baseU,baseUdash,baseUdashdash,gamma,Tb,B,C,shoot1,c,beta,tol)

    % Initialise vectors
    
    vec=[]; eigvec=[];

    % Loop through different khat values
    
    for shoot1=eigval-tol:tol/10:eigval+tol
    
        % Far field boudary condition 
        
        a1 = [exp(-((shoot1^2+beta^2)^0.5)*b),...
            -((shoot1^2+beta^2)^0.5)*exp(-((shoot1^2+beta^2)^0.5)*b)];
        %a2 = [exp(-khat*A^2/(3*a^3)), khat*A^2/(a^4)*exp(-khat*A^2/(3*a^3))];
        
        % Runge kutta inwards
        
        [eta, F1] = RK(a,b,deltaeta,a1,rayleigh,baseT,baseTdash,baseU,...
        baseUdash,baseUdashdash,gamma,Tb,shoot1,c,beta);
    
        % Boundary condition constraints
        
        H1=F1(2,1) + (9/a^4 + 2/a + B*((shoot1^2+beta^2)^0.5)/(a^(3-sqrt(7))) ...
            + 2*C*((shoot1^2+beta^2)^0.5))*F1(1,1);
        
        H2=F1(2,end) + ((shoot1^2+beta^2)^0.5)*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[H1,vec]; eigvec=[shoot1,eigvec];
    
    end
    
    p=F1;
    
% Calculate the crossing points
    
zerIdx=[];
for i=1:length(vec)-1
    if ((vec(i)>0 && vec(i+1)<0) || (vec(i)<0 && vec(i+1)>0))
    zerIdx(end+1)=i; % save index of zero-crossing
    end
end

eigs=eigvec(zerIdx);
vecs=vec(zerIdx);
eigval=eigs(1);
vec=vecs(1);