function x = newtonSolver(x, u, b, deltaT, a0, b0, Func, jacobF)

condRes = 1;
condDelta = 1;
absTol = 10^-12; % absolute tolerance
F = a0*x-deltaT*b0*Func(x,u)+b;
thrRes = absTol + 10^-6*norm(F); % threshold restriction
thrDelta = 0;

% Newton's method calculation
while condRes > thrRes || condDelta > thrDelta % convergence condition
    J = a0*eye(length(x))-deltaT*b0*jacobF(u, x); % Jacobian matrix calculation
    F = a0*x-deltaT*b0*Func(x,u)+b; % rhs calculation
    deltaX = J\-F; % solving J*deltaX = -F
    x = x + deltaX; % updating x at each step until convergence
    condDelta = norm(deltaX);
    thrDelta = absTol + 10^-6*norm(deltaX);
    condRes = norm(F);
end