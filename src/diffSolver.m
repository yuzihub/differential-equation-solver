function x = diffSolver(x0, deltaT, T, k, alpha, beta, Func, jacobF, u)

x = zeros(length(x0), T/deltaT+1); % declaring x
x(:,1) = x0; % setting the initial value of x

% implementation of the multi-step equation
for l = 2:T/deltaT+1
    
    b1 = zeros(size(x0));
    b2 = zeros(size(x0));
    if l <= k % first k step calculation
        n = 1;
        alpha1 = [1 -1];
        beta1 = [0.5 0.5];
        for j = 2:n+1 % redundant for loop kept in place for consistency
            b1 = b1 + alpha1(j)*x(:,l+1-j); % calculation of beta 1
            if size(u,2)>1 % calculation of beta 2
                b2 = b2 + beta1(j)*Func(x(:,l+1-j),u(l+1-j,:));
            else
                b2 = b2 + beta1(j)*Func(x(:,l+1-j),u(l+1-j));
            end
        end
    else % steps>k calculation
        n = k; 
        for j = 2:n+1
            b1 = b1 + alpha(j)*x(:,l+1-j);
            if size(u,2)>1
                b2 = b2 + beta(j)*Func(x(:,l+1-j),u(l+1-j,:));
            else
                b2 = b2 + beta(j)*Func(x(:,l+1-j),u(l+1-j));
            end
        end
    end
    

    b = b1 + -deltaT*b2;

    if size(u,2)>1
        x(:,l) = newtonSolver(x(:,l-1), u(l,:), b, deltaT, alpha(1), beta(1), Func, jacobF);
    else
        x(:,l) = newtonSolver(x(:,l-1), u(l), b, deltaT, alpha(1), beta(1), Func, jacobF);
    end
    
end