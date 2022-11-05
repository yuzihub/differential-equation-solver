function [alpha, beta] = coeffCalc(k)

array = zeros((k+1)*2-1, (k+1)*2);
for i = 1:(k+1)*2-1
    for j = 1:k+1
        array(i,j) = (k+1-j)^(i-1);
        array(i,j+k+1) = -(i-1)*(k+1-j)^(i-2);
    end
end

array(1,end) = 0;
array1 = array(:,2:end);
rhs = array(:,1);
coeffs = array1\-rhs;
alpha = cat(1, 1, coeffs(1:k));
beta = coeffs(k+1:end);