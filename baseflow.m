function [eta,baseT,baseTdash,baseU,baseUdash,baseUdashdash] = ...
    baseflow(C,Pr,D,etab,deltaeta,a,b)

dydx=@(eta,z)[z(2);z(3);(-eta*z(3)+ (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(3))/(((1+C)*sqrt(z(4)))/(z(4)+C)); ...
    z(5);(-eta*z(5) + (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5))/((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))];

BC=@(za,zb)[za(1) - D/(etab^(3/Pr)) ; zb(2) ; za(2) + (3/Pr)*D/etab^((3/Pr)-1); za(4)- (9*(1+C)^2)/(Pr^2*etab^4); zb(4)-1];
    
zint=@(x)[0 ; 1; 0 ; 1 ; 0];
    
solint=bvpinit(a:deltaeta:b,zint);
    
S=bvp4c(dydx,BC,solint);
    
eta=S.x; baseT=S.y(4,:); baseTdash=S.y(5,:);
 
baseU=S.y(2,:); baseUdash=S.y(3,:); baseUdashdash=S.y(3,:);

for i=2:length(baseUdash)-1
    baseUdashdash(i) = (baseUdash(i+1)-baseUdash(i-1))/(eta(2)-eta(1));
end

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

    