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



function [eta,baseU,baseUdash] = ...
    baseflow(C,Pr,deltaeta,a,b)

% Orig
%dydx = @(eta,z)[  z(2)  ;  z(3)  ;  ( -eta*z(3) + ...
%    (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2)) ...
%    *z(3) )/(((1+C)*sqrt(z(4)))/(z(4)+C))  ; ...
%    z(5)  ;  ( -eta*z(5) + ...
%    (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5)) ...
%    /((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))  ];

gamma=1.4; Tb=1;
c1=0.5*(gamma-1)*(Tb+1);

dydx=@(eta,z)[   z(2)   ;   z(3)   ;   ...
    ( -eta*z(3)*(C-c1*z(2)+1) )/( 2*((1-c1*z(2))^0.5) ) ...
    + ( c1*z(3)^2 )/( 2*(1-c1*z(2)) ) ...
    - ( 2*c1*z(3) )/( 2*(C-c1*z(2)+1) )   ];

% Boundary conditions

BC = @(za,zb)[  za(1) - 24/(2*c1*a^3)  ;  zb(1)  ; ...
    zb(2)  ];

% Initial conditions
    
zint = @(z)[ 0 ;  0  ;   0];

% Initialise solution

solint=bvpinit(a:deltaeta:b,zint);

% Solve with bvp4c

S=bvp4c(dydx,BC,solint);


% Extract solutions for base flow

eta=S.x; baseU=S.y(2,:); baseUdash=S.y(3,:);

% Interpolate for right grid size

baseU = interp1(eta,baseU,a:deltaeta:b,'spline');
baseUdash = interp1(eta,baseUdash,a:deltaeta:b,'spline');
eta=a:deltaeta:b;

% figure('position', [0,0,800,800]); 
% plot(eta,baseT,'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
% xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% 


% 
% figure('position', [0,0,800,800]); 
% plot(eta,baseU,'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Temp. adj. func, $G$','Interpreter', 'LaTex','Fontsize',40)
% xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on

    