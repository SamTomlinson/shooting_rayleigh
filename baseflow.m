function [x,baseT,baseTdash,baseU]= baseflow(C,Pr,D,eta)
    
dydx=@(x,z)[z(2);z(3);(-x*z(3)+ (((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(3))/(((1+C)*sqrt(z(4)))/(z(4)+C)); ...
    z(5);(-x*z(5) + (Pr^-1)*(((C+1)*(C-z(4)*z(5)))/(2*sqrt(z(4))*(C+z(4))^2))*z(5))/((Pr^-1)*((1+C)*sqrt(z(4)))/(z(4)+C))];

BC=@(za,zb)[za(1) - D/(eta^(3/Pr)) ; zb(2) ; za(2) + (3/Pr)*D/eta^((3/Pr)-1); za(4)- (9*(1+C)^2)/(Pr^2*eta^4); zb(4)-1];
    
zint=@(x)[0 ; 1; 0 ; 1 ; 0];
    
solint=bvpinit(linspace(1,3,1001),zint);
    
S=bvp4c(dydx,BC,solint);
    
x=S.x; baseT=S.y(4,:); baseTdash=S.y(5,:);
 
baseU=S.y(2,:);
    