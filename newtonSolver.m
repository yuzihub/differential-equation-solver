function x = newtonSolver(x, u, b, deltaT, a0, b0, Func, jacobF)

condRes = 1;
condDelta = 1;
absTol = 10^-12;
F = a0*x-deltaT*b0*Func(x,u)+b;
thrRes = absTol + 10^-6*norm(F);
thrDelta = 0;


while condRes > thrRes || condDelta > thrDelta
    J = a0*eye(length(x))-deltaT*b0*jacobF(u, x);
    F = a0*x-deltaT*b0*Func(x,u)+b;
    deltaX = J\-F;
    x = x + deltaX;
    condDelta = norm(deltaX);
    thrDelta = absTol + 10^-6*norm(deltaX);
    condRes = norm(F);
end