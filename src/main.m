clear all
close all

Func = @(x, u) [x(2); u*(1-x(1)^2)*x(2)-x(1)]; % the differential function as a system of first order ODEs
jacobF = @(u, x) [[0, 1]; [-2*u*x(2)*x(1)-1, u*(1-x(1)^2)]]; % the Jacobian matrix of the system 



k = 1; % number of steps used in a k-step multisteo method
x0 = [2; 0]; % initial values x(0) and x'(0)
deltaT = 0.01; % timestep
const = 1;
T = max(20,10*const); % end of time interval
interval = 0:deltaT:T; % time interval
betaArr = [0 1 1 0 0.5 0.5]; % beta coefficients for forward-Euler, backward-Euler, and trapezoidal rule
h =[];

for i = 1:3
    alpha = [1 -1]; % alpha0 and alpha1 coefficients
    beta = betaArr(i+(i-1)*1:i+1+(i-1)*1); % beta coefficient for each method
    
    u = const*ones(length(interval),1);
    
    x = diffSolver(x0, deltaT, T, k, alpha, beta, Func, jacobF, u); % solved diff. eqn.
    

    h = [h plot(interval, x)];
    hold on
end


% plotting the solutions for three different methods
title('Van der Pol Equation')
xlabel('Time (t)')
ylabel('Solution')
grid on
grid minor

set(h(1,1), 'Color', [255 150 0]/255);
set(h(2,1), 'Color', 'r');
set(h(1,2), 'Color', [0 0 178]/255);
set(h(2,2), 'Color', 'b');
set(h(1,3), 'Color', [0 150 0]/255);
set(h(2,3), 'Color', 'g');


leg1 = legend(h(:,1),'x_2','x_1');
ah1 = axes('position', get(gca,'position'), 'visible', 'off');
leg2 = legend(ah1, h(:,2),'x_2','x_1');
ah2 = axes('position', get(gca,'position'), 'visible', 'off');
leg3 = legend(ah2, h(:,3),'x_2','x_1');
title(leg1, 'Forward Euler');
title(leg2, 'Backward Euler');
title(leg3, 'Trapezoidal Rule');


