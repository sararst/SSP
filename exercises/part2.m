sigma_squared = 1;
A1 = sqrt(2);
phi1 = 0;
f1 = 1/8;
n = 32;
num_realizations = 100;

% Initialize arrays to store estimates
sigma_squared_estimates = zeros(1, num_realizations);
A1_estimates = zeros(1, num_realizations);
phi1_estimates = zeros(1, num_realizations);
f1_estimates = zeros(1, num_realizations);

% Generate 100 realizations
for i = 1:num_realizations
    % Generate sinusoid in white Gaussian noise
    k = 0:n-1;
    sk = A1 * cos(2 * pi * f1 * k + phi1);
    vk = sqrt(sigma_squared) * randn(1, n);
    yk = sk + vk;

    % Perform the DFT
    Yf = fft(yk);

    % Find the index with maximum magnitude in the frequency domain
    [~, l] = max(abs(Yf));
    
    % ML estimates
    f1_estimate = (l - 1) / n;  % Adjust for MATLAB indexing
    A1_estimate = 2 / n * abs(Yf(l));
    phi1_estimate = angle(Yf(l));
    sigma_squared_estimate = 1/n * sum((yk - A1_estimate * cos(2 * pi * f1_estimate * k + phi1_estimate)).^2);

    % Store the estimates
    f1_estimates(i) = f1_estimate;
    A1_estimates(i) = A1_estimate;
    phi1_estimates(i) = phi1_estimate;
    sigma_squared_estimates(i) = sigma_squared_estimate;
end

% Display the estimates
disp('ML Parameter Estimates:');
disp(['f1: Mean = ', num2str(mean(f1_estimates)), ', Std Dev = ', num2str(std(f1_estimates))]);
disp(['A1: Mean = ', num2str(mean(A1_estimates)), ', Std Dev = ', num2str(std(A1_estimates))]);
disp(['phi1: Mean = ', num2str(mean(phi1_estimates)), ', Std Dev = ', num2str(std(phi1_estimates))]);
disp(['sigma_squared: Mean = ', num2str(mean(sigma_squared_estimates)), ', Std Dev = ', num2str(std(sigma_squared_estimates))]);

% Display Cramer-Rao bounds
disp('Cramer-Rao Lower Bounds:');
disp(['CRB_sigma_squared: ', num2str(2 * sigma_squared^2 / n)]);
disp(['CRB_A1: ', num2str(2 * sigma_squared / n)]);
disp(['CRB_phi1: ', num2str(8 * sigma_squared / (n * A1^2))]);
disp(['CRB_f1: ', num2str(6 * sigma_squared / (pi^2 * n^3 * A1^2))]);