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
% adjustment region), Tb (wall temp), c (neutral phase speed)
%
% khat - combined wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE
%



%                         Output eigenfuntion                          %



function [eta, p] = shooting_rayleigh3(rayleigh,deltaeta,a,b,beta) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; Tb=1;
    B=-1.81447; c=-0.993937;

     % Solve for the base flow 
    
    [~,baseU,baseUdash] ...
        = baseflow(C,Pr,deltaeta,a,b);
    
    baseU;
    
    % Begin time

    tic;
    
    % Initialise vectors
    
    vec=[]; eigvec=[];

    % Loop through different alpha values
    
    for shoot1=0.1:0.01:0.4
    
        % Far field boudary condition 
        w=-((shoot1^2+beta^2)^0.5);
        a1 = [exp(w*b),w*exp(w*b)];
        a2 = [exp(w*(b+deltaeta)),w*exp(w*(b+deltaeta))];
        a3 = [exp(w*(b+2*deltaeta)),w*exp(w*(b+2*deltaeta))];
        a4 = [exp(w*(b+3*deltaeta)),w*exp(w*(b+3*deltaeta))];
        
        % Runge kutta inwards
        
        [~, F1] = AM(a,b,deltaeta,a1,a2,a3,a4,rayleigh,baseU,...
        baseUdash,gamma,Tb,shoot1,c,beta);
    
        % Boundary condition constraints
        
        H1=F1(2,1) + (-9/a^4 + 2/a + B*((shoot1^2 + beta^2)^0.5)/(a^(3-sqrt(7))) ...
            + 2*C*((shoot1^2+beta^2)^0.5))*F1(1,1);
        
        H2=F1(2,end) + ((shoot1^2 + beta^2)^0.5)*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[H1,vec]; eigvec=[shoot1,eigvec];
    
    end
    
    vec;
    p=F1;
    
% Plot H vs eig
    
    figure('position', [0,0,800,800]); 
    plot(eigvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{\alpha})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Wave no., $\hat{\alpha}$','Interpreter', 'LaTex','Fontsize',40)
    %ylim([-1,1])
    grid on
    hold off;
    
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

% Improve accuracy 
diff=1;
tol=0.01;
while abs(diff>1e-16)
    eigvalold=eigval;
    [eigval,H1]=loop(eigval,a,b,deltaeta,rayleigh,...
        baseU,baseUdash,gamma,Tb,B,C,c,beta,tol);
    diff=abs(eigvalold-eigval);
    tol=tol/10;
end

% Calculation of eigenmodes 

 % Far field boudary condition 
w=-((eigval^2+beta^2)^0.5);
a1 = [exp(w*b),w*exp(w*b)];
a2 = [exp(w*(b+deltaeta)),w*exp(w*(b+deltaeta))];
a3 = [exp(w*(b+2*deltaeta)),w*exp(w*(b+2*deltaeta))];
a4 = [exp(w*(b+3*deltaeta)),w*exp(w*(b+3*deltaeta))];        
% Runge kutta inwards
[eta, F1] = AM(a,b,deltaeta,a1,a2,a3,a4,rayleigh,baseU,...
        baseUdash,gamma,Tb,eigval,c,beta);
        
H1=F1(2,1) + (-9/a^4 + 2/a + B*((shoot1^2+beta^2)^0.5)/(a^(3-sqrt(7))) ...
            + 2*C*((eigval^2+beta^2)^0.5))*F1(1,1)
H2=F1(2,end) + ((eigval^2+beta^2)^0.5)*F1(1,end)
p=F1;

% Plotting of eigenomdes (if running evvsk % out)

figure('position', [0,0,800,800]); 
plot(eta,p(1,length(eta)),'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Vel. in the temp. adj. region $p_0$','Interpreter',...
      'LaTex','Fontsize',40)
xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([0.1,b])
grid on
hold off;

 
end