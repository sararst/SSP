% point iii
n = 1000;

Y = pingstats('isl.stanford.edu', n, 'v');

% Gaussian distribution
estimatedG_m = mean(Y);
estimatedG_sigma = std(Y);
disp(['Maximum Likelihood Estimator (MLE) for Gaussian distribution: mean = ', num2str(estimatedG_m)]);
disp(['Maximum Likelihood Estimator (MLE) for Gaussian distribution: standard deviation = ', num2str(estimatedG_sigma)]);


% Rayleigh distribution
estimatedR_sigma = sqrt((sum(Y.^2))/(2*n));
disp(['Maximum Likelihood Estimator (MLE) for Rayleigh distribution: sigma  = ', num2str(estimatedR_sigma)]);


% Erlang distribution
estimatedEm_lambda_val = zeros(1, 3);
for m = 0:2
    estimatedEm_lambda = n*(m+1)/sum(Y);
    estimatedEm_lambda_val(m+1) = estimatedEm_lambda; 
    disp(['Maximum Likelihood Estimator (MLE) for Erlang distribution with m = ', num2str(m), ':  ', num2str(estimatedEm_lambda)]);
end


% Shifted exponential distribution
estimatedExp_alpha = min(Y);
estimatedExp_lambda = n/(sum(Y-estimatedExp_alpha));
disp(['Maximum Likelihood Estimator (MLE) for shifted exponential distribution: alpha = ', num2str(estimatedExp_alpha)]);
disp(['Maximum Likelihood Estimator (MLE) for shifted exponential distribution: lambda = ', num2str(estimatedExp_lambda)]);


% Shifted Rayleigh
% Define the negative log-likelihood function
negative_log_likelihoodSR = @(alphaSR) -(-sum(1./(Y - alphaSR)) + (2 * n * sum(Y - alphaSR)) / sum((Y - alphaSR).^2));
% Initial guess for alpha
initial_alphaSR = min(Y) / 2;
lb = 0;
ub = min(Y);
% Optimize using fminsearch
[optimal_alphaSR, max_likelihoodSR] = fmincon(negative_log_likelihoodSR, initial_alphaSR, [], [], [], [], lb, ub, []);
% Display the result
disp(['Optimal alpha Shifted Rayleigh: ', num2str(optimal_alphaSR)]);
% Set initial guess for the scale parameter
alphaSR = optimal_alphaSR;
estimatedSR_sigma = sqrt(sum((Y-alphaSR).^2)/(2*n));
disp(['MLE for Shifted Rayleigh distribution: sigmaSR = ', num2str(estimatedSR_sigma)]);


% point iv
% calculate the range for y
y_min = min(Y);
y_max = max(Y);
y_range = linspace(y_min-20, y_max+20, 1000);

% Gaussian distribution
pdf_gaussian = normpdf(y_range, estimatedG_m, estimatedG_sigma);

% Rayleigh distribution
pdf_rayleigh = raylpdf(y_range, estimatedR_sigma);

% Erlang distribution
pdf_erlang0 = exppdf(y_range, 1/estimatedEm_lambda_val(1));
pdf_erlang1 = estimatedEm_lambda_val(2).^(2)*y_range.*exp(-estimatedEm_lambda_val(2)*y_range);
pdf_erlang2 = (estimatedEm_lambda_val(3).^(3)/2)*(y_range.^2).*exp(-estimatedEm_lambda_val(3)*y_range);

% Shifted exponential distribution
pdf_exponential_shifted = estimatedExp_lambda*exp(-estimatedExp_lambda*(y_range-estimatedExp_alpha));

% Shidted Rayleigh distribution
pdfSR = ((y_range - alphaSR) / estimatedSR_sigma^2).*exp(-((y_range - alphaSR).^2) / (2*estimatedSR_sigma.^2));

% plot histogram
figure;
histogram(Y, 'Normalization', 'probability', 'EdgeColor', 'w');
hold on;

% superimpose graphs for marginal densities
plot(y_range, pdf_gaussian, 'LineWidth', 2.5, 'DisplayName', 'Gaussian');
plot(y_range, pdf_rayleigh, 'LineWidth', 2.5, 'DisplayName', 'Rayleigh');
plot(y_range, pdf_erlang0, 'LineWidth', 2.5, 'DisplayName', 'Erlang (m=0)');
plot(y_range, pdf_erlang1, 'LineWidth', 2.5, 'DisplayName', 'Erlang (m=1)');
plot(y_range, pdf_erlang2, 'LineWidth', 2.5, 'DisplayName', 'Erlang (m=2)');
plot(y_range, pdf_exponential_shifted, 'LineWidth', 2.5, 'DisplayName', 'Shifted Exponential');
plot(y_range, pdfSR , 'LineWidth', 2.5, 'DisplayName', 'Shifted Rayleigh');

% add labels and legend
xlabel('Roundtrip delay');
ylabel('Probability density')
title('Histogram and superimposed densities')
legend('Histogram', 'Gaussian', 'Rayleigh', 'Erlang m=0', 'Erlang m=1', 'Erlang m=2', 'Shifted Exponential', 'Shifted Rayleigh')
grid on;

hold off;


% point v
log_likelihood_gaussian = -n*log(2*pi)/2 - n*log(estimatedG_sigma.^2)/2 -(sum(Y-estimatedG_m).^2)/(2*estimatedG_sigma.^2);
log_likelihood_rayleigh = sum(log(Y/estimatedR_sigma^2) - Y.^2/(2*estimatedR_sigma^2));
log_likelihood_erlang0 = sum(log(estimatedEm_lambda_val(1)) - estimatedEm_lambda_val(1)*Y);
log_likelihood_erlang1 = sum(2*log(estimatedEm_lambda_val(2))+log(Y)-estimatedEm_lambda_val(2)*Y);
log_likelihood_erlang2 = sum(3*log(estimatedEm_lambda_val(3)) - log(2) + 2*log(Y) -estimatedEm_lambda_val(3)*Y);
log_likelihood_shiftedExp = sum(log(estimatedExp_lambda)-estimatedExp_lambda*(Y-estimatedExp_alpha));
log_likelihood_SR = sum(log(Y - alphaSR) -log(estimatedSR_sigma^2) - (Y-alphaSR).^2/(2*estimatedSR_sigma.^2));


% Display the log-likelihoods
disp(['Log-Likelihood for Gaussian distribution: ', num2str(log_likelihood_gaussian)]);
disp(['Log-Likelihood for Rayleigh distribution: ', num2str(log_likelihood_rayleigh)]);
disp(['Log-Likelihood for Erlang distribution (m=0): ', num2str(log_likelihood_erlang0)]);
disp(['Log-Likelihood for Erlang distribution (m=1): ', num2str(log_likelihood_erlang1)]);
disp(['Log-Likelihood for Erlang distribution (m=2): ', num2str(log_likelihood_erlang2)]);
disp(['Log-Likelihood for shifted Exponential distribution: ', num2str(log_likelihood_shiftedExp)]);
disp(['Log-Likelihood for shifted Rayleigh distribution: ', num2str(log_likelihood_SR)]);


% List of distribution names
distribution_names = {'Gaussian', 'Rayleigh', 'Erlang (m=0)', 'Erlang (m=1)', 'Erlang (m=2)', 'Shifted Exponential', 'Shifted Rayleigh'};

log_likelihoods = [log_likelihood_gaussian, log_likelihood_rayleigh, log_likelihood_erlang0, log_likelihood_erlang1, log_likelihood_erlang2, log_likelihood_shiftedExp, log_likelihood_SR];

% Determine the best distribution
[best_log_likelihood, best_distribution_index] = max(log_likelihoods);

best_distribution_name = distribution_names{best_distribution_index};

disp(['The best distribution is: ', best_distribution_name]);
