%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 evvsk                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Runs shooting_rayeligh3 for a variety of different wavenumbers the plots
% the dependance of the ev on k. Plots some of the different eigenmodes



%                                 Key                                 % 
%
% eta - grid points
%
% p - p0 and p0dash array 
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
% eigval - found eigenvalue for entered khat
%
% k - spanwise wavenumber
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE



%                           Mostly plotting                           %



% Flow parameters 

gamma=1.4; Pr=1; C=0.509; Tb=1;
B=-1.81447; c=-0.993937;
deltaeta=0.01; a=0.1; b=5;

% Base flow

[eta,baseU,baseUdash] = baseflow(C,Pr,deltaeta,a,b);

plot(eta,baseU)

%Initialise 

ev=[];

% Loop through different wavenumbers and calculate ev

for beta=0.1:0.1:0.8
    
    beta
    
    % Calculate ev
    
    [eta, p,eigval] = shooting_rayleigh3(@rayleigh,deltaeta,a,b,beta);
    
    eigval
    
    % Fill vector
    
    ev=[eigval, ev];

    % Extract ones to plot
    
    if beta==0.1
        p1=p(1,:);
        dp1=p(2,:);
    end
    if beta==0.3
        p2=p(1,:); 
        dp2=p(2,:);
    end
    if beta==0.5
        p3=p(1,:);
        dp3=p(2,:);
    end     
end

% Reverse ev order as come out the wrong way round

ev=ev(end:-1:1);

%% Normalisation of other flow variables

% T1=normalise(-baseTdashm.*v1./baseTm);
% T2=normalise(-baseTdashm.*v2./baseTm); 
% T3=normalise(-baseTdashm.*v3./baseTm);
% u1=-baseUdashm.*v1./baseTm+1; 
% u2=-baseUdashm.*v2./baseTm+1; 
% u3=-baseUdashm.*v3./baseTm+1;
% w1=normalise(-dv1./(0.5.*baseTm)); 
% w2=normalise(-dv2./(1.*baseTm)); 
% w3=normalise(-dv3./(2.*baseTm));
% p1=normalise(-dv1./(0.5.^2.*baseTm.^2));
% p2=normalise(-dv2./(1.^2.*baseTm.^2)); 
% p3=normalise(-dv3./(2.^2.*baseTm.^2));


%% Plotting

% Recalculate base eta and base flow variables with correct size
% baseeta=a:deltaeta/5:b;
% baseTm = interp1(baseeta,baseT,a:deltaeta:b,'spline');
% baseTdashm = interp1(baseeta,baseTdash,a:deltaeta:b,'spline');
% baseUm = interp1(baseeta,baseU,a:deltaeta:b,'spline');
% baseUdashm = interp1(baseeta,baseUdash,a:deltaeta:b,'spline');
% 
% figure('position', [0,0,800,800]); 
% plot(baseeta,baseT,'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
% xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% 
% figure('position', [0,0,800,800]); 
% plot(baseeta,baseU,'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Vel. in adj. region, $U_1$','Interpreter', 'LaTex','Fontsize',40)
% xlabel('Wall layer variable, $\zeta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on


figure('position', [0,0,800,800]); 
plot(0.1:0.1:0.8,ev); hold on;
set(gca,'Fontsize',20)
ylabel('Wavenumber, $\alpha$','Interpreter',...
        'LaTex','Fontsize',40)
xlabel('Wavenumber, $\beta$','Interpreter', 'LaTex','Fontsize',40)
xlim([0.1,0.8])
grid on
hold off;

figure('position', [0,0,800,800]); 
plot(eta,p1,'LineWidth',2); hold on; 
plot(eta,p2,'LineWidth',2);
plot(eta,p3,'LineWidth',2);
set(gca,'Fontsize',20)
l1=legend('$\beta=0.1$','$\beta=0.3$','$\beta=0.5$');
set(l1, 'Interpreter','LaTex','Fontsize',30);
ylabel('Pressure in the temp. adj. region $p_0$','Interpreter',...
        'LaTex','Fontsize',40)
xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([a,b])

grid on
hold off;
%     
% figure('position', [0,0,800,800]); 
% plot(eta,T1,'LineWidth',2); hold on; 
% plot(eta,T2,'LineWidth',2); 
% plot(eta,T3,'LineWidth',2); 
% set(gca,'Fontsize',20)
% l1=legend('$k=0.5$','$k=1$','$k=2$');
% set(l1, 'Interpreter','LaTex','Fontsize',30);
% ylabel('Temp. in the temp. adj. region $T_0$','Interpreter',...
%         'LaTex','Fontsize',40)
% xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% hold off;
% 
% figure('position', [0,0,800,800]); 
% plot(eta,u1,'LineWidth',2); hold on; 
% plot(eta,u2,'LineWidth',2); 
% plot(eta,u3,'LineWidth',2); 
% set(gca,'Fontsize',20)
% l1=legend('$k=0.5$','$k=1$','$k=2$');
% set(l1, 'Interpreter','LaTex','Fontsize',30);
% ylabel('Vel. in the temp. adj. region $u_0$','Interpreter',...
%         'LaTex','Fontsize',40)
% xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% hold off;
% 
% figure('position', [0,0,800,800]); 
% plot(eta,w1,'LineWidth',2); hold on; 
% plot(eta,w2,'LineWidth',2); 
% plot(eta,w3,'LineWidth',2); 
% set(gca,'Fontsize',20)
% l1=legend('$k=0.5$','$k=1$','$k=2$');
% set(l1, 'Interpreter','LaTex','Fontsize',30);
% ylabel('Vel. in the temp. adj. region $w_0$','Interpreter',...
%         'LaTex','Fontsize',40)
% xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% hold off;
% 
% figure('position', [0,0,800,800]); 
% plot(eta,p1,'LineWidth',2); hold on; 
% plot(eta,p2,'LineWidth',2); 
% plot(eta,p3,'LineWidth',2); 
% set(gca,'Fontsize',20)
% l1=legend('$k=0.5$','$k=1$','$k=2$');
% set(l1, 'Interpreter','LaTex','Fontsize',30);
% ylabel('Pres. in the temp. adj. region $p_0$','Interpreter',...
%         'LaTex','Fontsize',40)
% xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% hold off;