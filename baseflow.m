%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               baseflow                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%                           Code description                          %



% Solves boundary value problem for the base flow




%                                 Key                                 % 
%
% eta - grid points
%
% deltaeta - step size
%
% zint - initial guess for solution profiles
%
% a,b - two ends of the domain
%
% flow parameters - gamma (specific heat), Pr (prandtl), C (sutherlands
% constant), D (fitting parameter), etab (matching point for edge of
% adjustment region)
% 
% base flow - baseT, baseTdash, base U base flow vectors and derivatives
% intbaseT integral for Q term in DE
%



%                         Output base flow solve                        %



function [eta,baseT,baseTdash,baseU,baseUdash,baseUdashdash] = ...
    baseflow(C,Pr,D,etab,deltaeta,a,b)

% Orig
dydx = @(eta,z)[  z(2)  ;  z(3)  ;  ( -eta*z(3) + ...
    (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2)) ...
    *z(3) )/(((1+C)*sqrt(z(4)))/(z(4)+C))  ; ...
    z(5)  ;  ( -eta*z(5) + ...
    (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5)) ...
    /((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))  ];

% New
% dydx = @(eta,z)[  z(2)  ;  z(3)  ;  ...
%     (-eta*z(3)*(z(4)+C))/(sqrt(z(4))*(1+C)) ...
%     - ((C-z(4))*z(5)^2)/(2*z(4)*(C+z(4)))  ; ...
%     z(5)  ;  (-Pr*eta*z(5)*(C+z(4)))/((1+C)*sqrt(z(4))) - ...
%     ((C-z(4))*z(3)*z(5))/(2*z(4)*(C+z(4)))];

% Boundary conditions

BC = @(za,zb)[  za(1) - D/(etab^(3/Pr))  ;  zb(2)  ; ...
    za(2) + (3/Pr)*D/etab^((3/Pr)-1)  ; ...
    za(4)- (9*(1+C)^2)/(Pr^2*etab^4)  ;  zb(4)-1  ];

% Initial conditions
    
zint = @(x)[  0  ;  1  ;  0  ;  1  ;  0 ];

% Initialise solution

solint=bvpinit(a:deltaeta*10:b,zint);
    
% Solve with bvp4c

S=bvp4c(dydx,BC,solint);

% Extract solutions for base flow
    
eta=S.x; baseT=S.y(4,:); baseTdash=S.y(5,:);
baseU=S.y(2,:); baseUdash=S.y(3,:); baseUdashdash=S.y(3,:);
for i=2:length(baseUdash)-1
    baseUdashdash(i) = (baseUdash(i+1)-baseUdash(i-1))/(eta(2)-eta(1));
end
baseUdashdash(1)=(3/Pr)*((3/Pr)-1)*D/etab^((3/Pr)-2);
baseUdashdash(end)=0;

% Interpolate for right grid size

baseT = interp1(eta,baseT,a:deltaeta:b,'spline');
baseTdash = interp1(eta,baseTdash,a:deltaeta:b,'spline');
baseU = interp1(eta,baseU,a:deltaeta:b,'spline');
baseUdash = interp1(eta,baseUdash,a:deltaeta:b,'spline');
baseUdashdash = interp1(eta,baseUdashdash,a:deltaeta:b,'spline');
eta=a:deltaeta:b;

figure('position', [0,0,800,800]); 
plot(eta,baseT,'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([a,b])
grid on

figure('position', [0,0,800,800]); 
plot(eta,baseU,'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Vel. in adj. region, $U_1$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\zeta$','Interpreter', 'LaTex','Fontsize',40)
xlim([a,b])
grid on

    