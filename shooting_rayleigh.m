%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           shooting_gotler                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Code description                          %

% This code implements the shooting method for solving 1D boundary value 
% problem. It uses the Runge-Kutta method of 4th order for solving the 
% ODE and the interval bisection method for finding the alpha parameter. 

% The meaning of input and output parameters: 
%
% x - grid of points where the solution have been found
%
% y - 2D array, values of function (solution) are in first row, values 
% of 1st derivative are in second row
%
% gotler - the function handle, fun contains system of solved 
% differential equations
%
% h - the step of the Runge-Kutta method (the step of the grid)            
% zero - the interval bisection method accuracy
%
% con - values of boundary conditions (2D vector)
%
% type - to specify type of the boundary condition in the particular
% point, it is the string consists of 2 char ('f' for function, 'd' for
% derivative), e.g 'fd' is meant that in a point the condition is
% related to function and in the b point to derivative
%
% init - 2D vector of initial alpha parameters, it is not required, the 
% implicit value is [-10 10]

% The ploting of the solution:
% The solution (both function and 1st derivative) of BVP ODE is shown 
% graphicaly after enumeration.

% Example run that works to an extent:
%
% [x y,baseT] = shooting_gotler(@gotler,0.0030,1e-6,1,3,[0 0],'ff');
%
% It is meant that will be solved the BVP ODE described in the function
% gotler, on the interval (1,3) with boundary conditions y(1) = 0 and 
% y(3) = 0 to a tolerance of 1e-6.

% Function for output of eigenfuntion

function [x, y, baseT] = shooting_gotler(gotler,h,zero,a,b,con,...
    type,init) 

    % Parameters and base flow should really be put into funtion 

    gamma=1.4; Pr=1; C=0.509;
    D=1; % Fitting parameter for base flow 
    eta=1; % Chosen matching point or left boundary 
    betag=1; Gstar=2; Q=1; sigma=1; kappa=1;
    
    % Solve for the base flow 
    
    [x,baseT,baseTdash,baseU]= baseflow(C,Pr,D,eta);

    tic; % Begin time
    
    % If my number of arguements is 8 then initial guesses gave been 
    % specified if not take these to be -1 and 1.
    
    if nargin == 8
        shoot1 = init(1); shoot2 = init(2);
    else
        shoot1 = -20; shoot2 = 10;
    end
    
    % Sets up boundary condition vectors, with the first entries being
    % the know dirichlet conditions and the second the two shoots
    
    if (type(1)=='f')
        a1 = [con(1) shoot1];
        a2 = [con(1) shoot2];
    else
        a1 = [shoot1 con(1)];
        a2 = [shoot2 con(1)];
    end  
    
    % Now iterate solution outwards using Rk method 
    
    [x, F1] = RK(a,b,h,a1,gotler,baseT,baseTdash,baseU,kappa,betag,...
        Gstar,Q,sigma); 
    [x, F2] = RK(a,b,h,a2,gotler,baseT,baseTdash,baseU,kappa,betag,...
        Gstar,Q,sigma);         
    
    if (type(2)=='f')
        F1 = F1(1,end) - con(2);
        F2 = F2(1,end) - con(2);
        r = 1;
    else
        F1 = F1(2,end) - con(2); 
        F2 = F2(2,end) - con(2);
        r = 2;
    end    
    
    % Identify if as root is possible by checking for sign change
    
    if (F1*F2 > 0) 
        error('The root of F function does not exist')
    end
    
    % Set one shoot for iteration 
    
    F3 = F1;
    
    % Iteration to home in on axis crossing 
    
    while (abs(F3) > zero) 
        
        % Bring one shoot in half the distance between the teo
        
        shoot3 = (shoot1 + shoot2)/2;
        
        % Renforce conditions and rerun RK solver on loop adjusting 
        % to compensate for average overshooting root
        
        if (type(1)=='f')
           a3 = [con(1) shoot3];            
        else
           a3 = [shoot3 con(1)];            
        end           
        
        [x, F3] = RK(a,b,h,a3,gotler,baseT,baseTdash,baseU,kappa,...
            betag,Gstar,Q,sigma);
        
        y = F3; F3 = F3(r,end) - con(2); 
        if (F1*F3 < 0)
            shoot2 = shoot3; F2 = F3;            
        elseif (F1*F2 < 0)
            shoot1 = shoot3; F1 = F3;
        else
            error('Something has gone horribly wrong, probs NANS');           
        end
    end           
    
    % Plotting of solutions 
    
    figure('position', [0,0,800,800]); 
    plot(x,y(1,:),'k-','LineWidth',2); hold on; 
    plot(x,y(2,:),'r-','LineWidth',2); 
    set(gca,'Fontsize',20)
    l1=legend('$v_0(\eta)$','$v_{0\eta}(\eta)$');
    set(l1, 'Interpreter','LaTex','Fontsize',30);
    ylabel('Vel. in the temp. adj. region $v_0$','Interpreter', 'LaTex','Fontsize',40)
    xlabel('D.H. variable, $\eta$','Interpreter', 'LaTex','Fontsize',40)
    xlim([1,3])
    grid on
    hold off;
    toc
    
end