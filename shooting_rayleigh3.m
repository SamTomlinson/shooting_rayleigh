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



function [eta, p,eigval] = shooting_rayleigh3(rayleigh,deltaeta,a,b,beta) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509; Tb=1;
    B=-1.81; c=-0.99;

     % Solve for the base flow 
    
    [~,baseU,baseUdash] ...
        = baseflow(C,Pr,deltaeta,a,b+4*deltaeta);
%     [~,baseU,baseUdash] ...
%         = baseflow(C,Pr,deltaeta,a,b);
    
    % Begin time

    tic;
    
    % Initialise vectors
    
    vec=[]; eigvec=[];

    % Loop through different alpha values
    
    for shoot1=0.1:0.1:1
    
        % Far field boudary condition extra for multistep methods
        w=-((shoot1^2+beta^2)^0.5);
        a1 = [exp(w*b),w*exp(w*b)];
        a2 = [exp(w*(b+deltaeta)),w*exp(w*(b+deltaeta))];
        a3 = [exp(w*(b+2*deltaeta)),w*exp(w*(b+2*deltaeta))];
        a4 = [exp(w*(b+3*deltaeta)),w*exp(w*(b+3*deltaeta))];
        
        % Iterate inwards
        
       [~, F1] = AM(a,b,deltaeta,a1,a2,a3,a4,rayleigh,baseU,...
        baseUdash,gamma,Tb,shoot1,c,beta);
%          [~, F1] = RK(a,b,deltaeta,a1,rayleigh,baseU,...
%          baseUdash,gamma,Tb,shoot1,c,beta);
    
        % Boundary condition constraints
        
        H1=F1(2,1) + ( (-9*((shoot1^2 + beta^2)^0.5))/a^4 + 2/a ...
            + (B*((shoot1^2 + beta^2)^0.5))/(a^(3-sqrt(7))) ...
            + 2*C*((shoot1^2+beta^2)^0.5) )*F1(1,1);
        
%         Hposs=exp(-(((shoot1^2 + beta^2)^0.5)*((B*a^(1+sqrt(7)))/(sqrt(7)-2) ...
%            + 2*C*a^4 + 3))/(a^3) - 2*log(a))
        
        H2=F1(2,end) + (shoot1^2 + beta^2)*F1(1,end);
   
        % Vector of H error and ks
        
        vec=[vec,H1]; eigvec=[eigvec,shoot1];
    
    end
    
    eigvec;
    vec;
    
% Plot H vs eig
    
    figure('position', [0,0,800,800]); 
    plot(eigvec,vec,'k-','LineWidth',2); 
    set(gca,'Fontsize',20)
    ylabel('Near field error, $H(\hat{\alpha})$','Interpreter',...
        'LaTex','Fontsize',40)
    xlabel('Wave no., $\hat{\alpha}$','Interpreter', 'LaTex','Fontsize',40)
    xlim([0,1])
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
eigval=eigs(1)
vec=vecs(1);

%Improve accuracy 
diff=1;
tol=deltaeta*10;
k=1;
while abs(k<16)
    k
    eigvalold=eigval
    [eigval,H1]=loop(eigvalold,a,b,deltaeta,rayleigh,...
        baseU,baseUdash,gamma,Tb,B,C,c,beta,tol);
    diff=abs(eigvalold-eigval);
    tol=tol/10;
    k=k+1;
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
    
% [eta, F1] = RK(a,b,deltaeta,a1,rayleigh,baseU,...
%         baseUdash,gamma,Tb,eigval,c,beta);
        
H1=F1(2,1) + (-9*((eigval^2 + beta^2)^0.5)/a^4 + 2/a ...
            + B*((eigval^2 + beta^2)^0.5)/(a^(3-sqrt(7))) ...
            + 2*C*((eigval^2+beta^2)^0.5))*F1(1,1);
H2=F1(2,end) + ((eigval^2+beta^2))*F1(1,end);
p=F1;

% Plotting of eigenomdes (if running evvsk % out)

% figure('position', [0,0,800,800]); 
% plot(eta,p(1:length(eta)),'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Pres. in the temp. adj. region $p_0$','Interpreter',...
%       'LaTex','Fontsize',40)
% xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% %ylim([0,1])
% grid on
% hold off;

vq1 = interp1(eta,p(1:length(eta)),a:10*deltaeta:b);
eta2=a:10*deltaeta:b;
figure()
plot(eta2,vq1)

end