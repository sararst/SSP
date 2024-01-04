sigma_squared = 1;
A1 = sqrt(2);
phi1 = 0;
f1 = 1/8;
n = 32;
num_realizations = 100;

% Initialize arrays to store estimates
sigma_squared_estimates = zeros(1, num_realizations);
A1_estimates = zeros(1, num_realizations);
f1_estimates = zeros(1, num_realizations);
realizations = zeros(n, num_realizations);

% Generate 100 realizations
for i = 1:num_realizations
    % Generate sinusoid in white Gaussian noise
    k = 0:n-1;
    sk = A1 * cos(2 * pi * f1 * k + phi1);
    vk = sqrt(sigma_squared) * randn(1, n);
    yk = sk + vk;
    % store realizations for point f
    realizations(:, i) = yk;
end

for i = 1:num_realizations
    % calculate sample autocorrelations
    r0 = mean(realizations(:, i).^2);
    r1 = mean(realizations(1:end-1, i) .* realizations(2:end, i));
    r2 = mean(realizations(1:end-2, i) .* realizations(3:end, i));

    % if r1 is different than 0
    if r1 ~= 0
        X = (r2+sqrt(r2.^2 + 8 * r1.^2))/(4 * r1);
        f1_estimates(i) = acos(X) / (2*pi);
        A1_estimates(i) = sqrt(2 * r1) / X;
        sigma_squared_estimates = r0 - (A1_estimates(i).^2)/2;
    else
        X = 0;
        f1_estimates(i) = 1/4;
        A1_estimates(i) = sqrt(-2 * r2);
        sigma_squared_estimates(i) = r0 - (A1_estimates(i).^2)/2;
    end
end

% Compute mean and variance for each estimate
f1_est_mean = mean(f1_estimates);
f1_est_var = var(f1_estimates);
A1_est_mean = mean(A1_estimates);
A1_est_var = var(A1_estimates);
sigma_sq_est_mean = mean(sigma_squared_estimates);
sigma_sq_est_var = var(sigma_squared_estimates);
disp('Covariance Matching Parameter Estimates:');
disp(['f1: Mean = ', num2str(f1_est_mean), ', Variance = ', num2str(f1_est_var)]);
disp(['A1: Mean = ', num2str(A1_est_mean), ', Variance = ', num2str(A1_est_var)]);
disp(['sigma_squared: Mean = ', num2str(sigma_sq_est_mean), ', Variance = ', num2str(sigma_sq_est_var)]);

% Cramer-Rao bounds
CRB_sigma_sq = 2 * sigma_squared^2 / n;
CRB_A1 = 2 * sigma_squared / n;
CRB_f1 = 6 * sigma_squared / (pi^2 * n^3 * A1^2);
disp('Cramer-Rao Lower Bounds:');
disp(['CRB_sigma_squared: ', num2str(CRB_sigma_sq)]);
disp(['CRB_A1: ', num2str(CRB_A1)]);
disp(['CRB_f1: ', num2str(CRB_f1)]);