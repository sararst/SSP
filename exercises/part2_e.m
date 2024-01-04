% point e
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

    % Compute the DFT
    Yl = fft(yk);

    % hat{l} = argmax|Y_l|
    [~, l] = max(abs(Yl));
    
    % ML estimates
    f1_estimate = (l - 1) / n;  % Adjust for MATLAB indexing
    A1_estimate = 2 / n * abs(Yl(l));
    phi1_estimate = angle(Yl(l));
    sigma_squared_estimate = 1/n * sum((yk - A1_estimate * cos(2 * pi * f1_estimate * k + phi1_estimate)).^2);

    % Store the estimates
    f1_estimates(i) = f1_estimate;
    A1_estimates(i) = A1_estimate;
    phi1_estimates(i) = phi1_estimate;
    sigma_squared_estimates(i) = sigma_squared_estimate;
end

% Compute mean and variance for each estimate
f1_est_mean = mean(f1_estimates);
f1_est_var = var(f1_estimates);
A1_est_mean = mean(A1_estimates);
A1_est_var = var(A1_estimates);
phi1_est_mean = mean(phi1_estimates);
phi1_est_var = var(phi1_estimates);
sigma_sq_est_mean = mean(sigma_squared_estimates);
sigma_sq_est_var = var(sigma_squared_estimates);
disp('ML Parameter Estimates:');
disp(['f1: Mean = ', num2str(f1_est_mean), ', Variance = ', num2str(f1_est_var)]);
disp(['A1: Mean = ', num2str(A1_est_mean), ', Variance = ', num2str(A1_est_var)]);
disp(['phi1: Mean = ', num2str(phi1_est_mean), ', Variance = ', num2str(phi1_est_var)]);
disp(['sigma_squared: Mean = ', num2str(sigma_sq_est_mean), ', Variance = ', num2str(sigma_sq_est_var)]);

% Cramer-Rao bounds
CRB_sigma_sq = 2 * sigma_squared^2 / n;
CRB_A1 = 2 * sigma_squared / n;
CRB_phi1 = 8 * sigma_squared / (n * A1^2);
CRB_f1 = 6 * sigma_squared / (pi^2 * n^3 * A1^2);
disp('Cramer-Rao Lower Bounds:');
disp(['CRB_sigma_squared: ', num2str(CRB_sigma_sq)]);
disp(['CRB_A1: ', num2str(CRB_A1)]);
disp(['CRB_phi1: ', num2str(CRB_phi1)]);
disp(['CRB_f1: ', num2str(CRB_f1)]);

