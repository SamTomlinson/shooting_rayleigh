function vecq = rayleigh(~,q,baseT,~,baseU,baseUdash,baseUdashdash,...
    c,gamma,Tb,khat)

    vecq(1) = q(2);
    
    % The %'d out option is the full equation
    
    %vecq(2) = 2*baseUdash*q(2)/(baseU-c) ...
    %    + ktilde^2*baseT*(baseT-alpha^2*M^2*(baseU-c)^2/ktilde^2)*q(1);
    
    % This option is the hypersonic limit
    
    vecq(2) = 2*(0.5*(gamma-1)*(Tb+1)*baseUdashdash*q(2))/...
        (0.5*(gamma-1)*(Tb+1)*baseUdash-c) ...
        + khat^2*(1-(0.5*(gamma-1)*(Tb+1)*baseUdash)^2)*q(1);
    
    