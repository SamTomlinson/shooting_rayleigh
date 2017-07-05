function vecq = gotler(~,q,baseT,baseTdash,baseU,kappa,betag,~,~,sigma)
    vecq(1) = q(2);
    vecq(2) = q(2)*baseTdash./baseT ...
        + ((baseT^(-2)).*betag.*q(1) + ...
        (0.5*kappa*(baseU^2)-trapz(baseT))*(baseTdash*betag/sigma))*q(1);
    
    