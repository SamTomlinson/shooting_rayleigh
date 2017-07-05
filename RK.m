function [x, y] = RK(a,b,h,con,rayleigh,baseT,baseTdash,baseU,baseUdash,...
    baseUdashdash,alpha,ktilde,M,c,gamma,Tb)
    x = a:h:b; 
    n = length(x);
    y = zeros(length(con),n);
    y(:,1) = con;
    
    for k = 1:(n-1) 
        k1 = rayleigh(x(k),y(:,k),baseT(k),baseTdash(k),baseU(k),...
        baseUdash(k),baseUdashdash,alpha,ktilde,M,c,gamma,Tb);
        k2 = rayleigh(x(k)+0.5*h,y(:,k)+0.5*h*k1',baseT(k),...
            baseTdash(k),baseU(k),baseUdash(k),baseUdashdash,...
            alpha,ktilde,M,c,gamma,Tb);
        k3 = rayleigh(x(k)+0.5*h,y(:,k)+0.5*h*k2',baseT(k),...
            baseTdash(k),baseU(k),baseUdash(k),baseUdashdash,...
            alpha,ktilde,M,c,gamma,Tb);
        k4 = rayleigh(x(k)+h,y(:,k)+h*k3',baseT(k),baseTdash(k),...
            baseU(k),baseUdash(k),baseUdashdash,alpha,ktilde,M,c,...
            gamma,Tb);
        y(:,k+1)=y(:,k)+h*(k1'+2*k2'+2*k3'+k4')/6;
    end