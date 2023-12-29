function [sigma2_f, A, K] = levinson(correlation_sequence)
% Input:
% correlation sequence r_{0:20}
% Output:
% sigma2_f = sequence of prediction error variances σ^2_{f,0:20}
% A = prediction error filter coefficients of order 20 A_{20, 0:20}
% K = sequence of PARCORS K_{1:20}
% Initialization
 order = length(correlation_sequence) - 1;

    % Initialize outputs
    sigma2_f = zeros(1, order + 1);
    A = zeros(order + 1, order + 1);
    K = zeros(1, order);

    % Initialize the recursion
    sigma2_f(1) = correlation_sequence(1);

    for n = 1:order
        % Compute the PARCOR
        K(n) = -correlation_sequence(n + 1) / sigma2_f(n);

        % Update sigma2_f
        sigma2_f(n + 1) = sigma2_f(n) * (1 - K(n)^2);

        % Update A coefficients using Levinson recursion
        A(n + 1, n + 1) = K(n);
        A(n + 1, 1:n) = A(n, 1:n) + K(n) * flip(A(n, 1:n));
    end
end