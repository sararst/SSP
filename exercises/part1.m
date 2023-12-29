% point B

% generate one realization of yk with N = 256 samples
N = 256;
[y, ~] = sig(N);

% define the window
window_type = 'boxcar';

% define the range of N' values (called Np)
Np_values = [64, 128, 256, 512, 1024];

figure;

for i = 1:length(Np_values)
    Np = Np_values(i);
    win = window(window_type, Np);
    
    % apply the window to the signal
    % y_windowed = y(1:min(length(y), length(win))) .* win(1:min(length(y), length(win)));
    if Np < N
        y_windowed = y(1:Np) .* win(1:Np);
    else
        y_windowed = [y; zeros(Np - N, 1)] .* win;
    end
    subplot(2, 3, i);
    % compute the periodogram
    periodo(y_windowed, Np);
    title(['Periodogram with Np = ' num2str(Np)]);
end

sgtitle('Effect of N prime on periodogram');



% point C
N = length(y);

% Reverse the signal
y_reverse = y(N:-1:1);

% Convolve the signal with its reversed version
correlation_result = conv(y, y_reverse);

% Extract the correlation sequence r0:n
n_max = min(N - 1, length(correlation_result) - N + 1);
correlation_sequence = correlation_result(N:N + n_max) / N;

% Display the correlation sequence
figure;
stem(0:n_max, correlation_sequence);
title('Correlation Sequence');
xlabel('n');
ylabel('Correlation Value');
grid on;


% point D

prediction_order = 20;
% Call the levinson function
[sigma2_f, A, K] = levinson(correlation_sequence(1:prediction_order+1));


% autoregressive spectrum estimate
figure;
N = 1024; 
A_padded = [1; A(:, end); zeros(N - length(A), 1)];
AR_spectrum = abs(fft(A_padded, N)).^2;
freq_axis = (0:N - 1) / N;


% Plot the periodogram
subplot(2, 1, 1);
periodo(y, 1024);
title('Periodogram with N = 1024');

% Plot the autoregressive spectrum estimate
subplot(2, 1, 2);
plot(freq_axis, 10 * log10(AR_spectrum));
hold on;

title('Autoregressive Spectrum Estimate');
xlabel('Normalized Frequency');
ylabel('dB');
grid on;
% Adjust layout
sgtitle('Periodogram and Autoregressive Spectrum Estimate');


% Plot the values of A20,0:20
subplot(2, 1, 2);
stem(0:prediction_order, A(:,end));
title('Values of A20,0:20');
xlabel('Coefficient Index');
ylabel('Coefficient Value');
grid on;

% Plot the evolution of σ^2_f,0:20 and the PARCORS
figure;

subplot(2, 1, 1);
stem(0:prediction_order, sigma2_f);
title('Evolution of \sigma^2_f,0:20');
xlabel('Order');
ylabel('\sigma^2_f');
grid on;

subplot(2, 1, 2);
stem(1:prediction_order, K);
title('Evolution of PARCORS K1:20');
xlabel('Order');
ylabel('PARCORS');
grid on;

