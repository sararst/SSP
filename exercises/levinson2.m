function [sigma_square, A, K] = levinson2(r)

% Initialization
p = length(r) - 1;
sigma_square = zeros(p+1, 1);
A = zeros(p+1, p+1);
K = zeros(p, 1);

% Step 0
sigma_square(1) = r(1);
A(1, 1) = 1;

% Step 1
K(1) = -r(2) / sigma_square(1);
A(2, 1) = K(1);
A(2, 2) = 1 + K(1);

% Step 2
for i = 2:p
    % Compute sigma_square(i)
    sigma_square(i) = sigma_square(i-1) * (1 - K(i-1)^2);
    
    % Compute K(i)
    K(i) = -r(i+1);
    for j = 1:i-1
        K(i) = K(i) - A(i+1, j+1) * r(i-j+1);
    end
    K(i) = K(i) / sigma_square(i);
    
    % Compute A(i+1, :)
    A(i+1, 1) = K(i);
    for j = 1:i
        A(i+1, j+1) = A(i, j) + K(i) * A(i, i-j+2);
    end
    A(i+1, i+2) = 1;
end

% Compute sigma_square(p+1)
sigma_square(p+1) = sigma_square(p) * (1 - K(p)^2);
