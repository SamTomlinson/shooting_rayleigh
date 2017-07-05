function vecq = rayleigh(~,q,baseT,~,baseU,baseUdash,alpha,ktilde,M,c)
    vecq(1) = q(2);
    vecq(2) = 2*baseUdash*q(2)/(baseU-c) ...
        + ktilde^2*baseT*(baseT-alpha^2*M^2*(baseU-c)^2/ktilde^2)*q(1);
    
    