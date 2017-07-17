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



function [eta,baseU,baseUdash,baseUdashdash,solint] = ...
    baseflow(C,Pr,D,deltaeta,a,b)

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
    za(2) + 3*24./(2*c1*a^4)   ];

% Initial conditions
    
zint = @(x)[   10 ;  0  ;   0];

% Initialise solution

solint=bvpinit(a:deltaeta/10:b,zint);

% Solve with bvp4c

S=bvp4c(dydx,BC,solint);

% dydx=@(x,z)[z(2);z(3);(-x.*z(3).*(C-c1.*z(2)+1))./(2.*(1-c1.*z(2)).^0.5) ...
%     + (c1.*z(3))./(2.*(1-c1.*z(2))) ...
%     - (c1.*z(3))./(C-c1.*z(2)+1)];
% BC=@(za,zb)[za(1) - 24./(2.*c1.*2^3) ; zb(1) ; za(2)+ (3).*24./(c1.*etab.^(4))];
% zint=@(x)[0 ; 0; 0];
% solint=bvpinit(linspace(2,15,100),zint);
% S=bvp4c(dydx,BC,solint);


% Extract solutions for base flow
    
eta=S.x; 
%baseT=S.y(4,:); baseTdash=S.y(5,:);
baseU=S.y(2,:); baseUdash=S.y(3,:); baseUdashdash=S.y(3,:);
for i=2:length(baseUdash)-1
    baseUdashdash(i) = (baseUdash(i+1)-baseUdash(i-1))/(eta(2)-eta(1));
end
baseUdashdash(1)=-3*4*24/(c1*a^(5));
baseUdashdash(end)=0;

% Interpolate for right grid size

%baseT = interp1(eta,baseT,a:deltaeta/5:b,'spline');
%baseTdash = interp1(eta,baseTdash,a:deltaeta/5:b,'spline');

baseU = interp1(eta,baseU,a:deltaeta/5:b,'spline');
baseUdash = interp1(eta,baseUdash,a:deltaeta/5:b,'spline');
baseUdashdash = interp1(eta,baseUdashdash,a:deltaeta/5:b,'spline');
eta=a:deltaeta/5:b;

% figure('position', [0,0,800,800]); 
% plot(eta,baseT,'LineWidth',2); 
% set(gca,'Fontsize',20)
% ylabel('Temp. in adj. region, $T_1$','Interpreter', 'LaTex','Fontsize',40)
% xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
% xlim([a,b])
% grid on
% 
figure('position', [0,0,800,800]); 
plot(eta,c1*baseU,'LineWidth',2); 
set(gca,'Fontsize',20)
ylabel('Temp. adj. func, $G$','Interpreter', 'LaTex','Fontsize',40)
xlabel('Wall layer variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
xlim([a,b])
grid on

    