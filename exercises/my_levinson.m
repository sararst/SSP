function [sigma2_f, A, K] = my_levinson(correlation_sequence)
    % Input:
    % correlation sequence r_{0:20}
    % Output:
    % sigma2_f = sequence of prediction error variances sigma^2_{f,0:20}
    % A = prediction error filter coefficients of order 20 A_{20, 0:20}
    % K = sequence of PARCORS K_{1:20}

    n=20;
    K = zeros(1, n+1);
    sigma2_f = zeros(1, n);

    % Initialization
    A = 1;
    sigma2_f(1) = correlation_sequence(1);

    for i=1:n
        delta = correlation_sequence(i+1:-1:2)*A';
        K(i+1) = -delta/sigma2_f(i);
        A = [A, 0] + K(i+1)*[0, A(i:-1:1)];
        sigma2_f(i+1) = sigma2_f(i)*(1-K(i+1)^2);
    end
    K = K(2:n+1)';