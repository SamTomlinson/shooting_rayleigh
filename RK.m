function [x, y] = RK(a,b,h,con,gotler,baseT,baseTdash,baseU,kappa,betag,Gstar,Q,sigma)
    x = a:h:b; 
    n = length(x);
    y = zeros(length(con),n);
    y(:,1) = con;
    
    for k = 1:(n-1) 
        k1 = gotler(x(k),y(:,k),baseT(k),baseTdash(k),baseU(k),kappa,betag,...
            Gstar,Q,sigma);
        k2 = gotler(x(k)+0.5*h,y(:,k)+0.5*h*k1',baseT(k),...
            baseTdash(k),baseU(k),kappa,betag,Gstar,Q,sigma);
        k3 = gotler(x(k)+0.5*h,y(:,k)+0.5*h*k2',baseT(k),...
            baseTdash(k),baseU(k),kappa,betag,Gstar,Q,sigma);
        k4 = gotler(x(k)+h,y(:,k)+h*k3',baseT(k),baseTdash(k),...
            baseU(k),kappa,betag,Gstar,Q,sigma);
        y(:,k+1)=y(:,k)+h*(k1'+2*k2'+2*k3'+k4')/6;
    end