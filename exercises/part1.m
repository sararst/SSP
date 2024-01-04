% point B
% generate one realization of yk with N = 256 samples
N = 256;
[y, ~] = sig(N);

% define the window
window_type = 'boxcar';
w = window(window_type, N);
y_windowed = y.*w;
figure;
subplot(2, 3, 1)
plot(y_windowed)
xlabel('Signal')
ylabel('y')
grid on

% define the range of N' values (called Np)
Np_values = [64, 128, 256, 512, 1024];

for i = 1:length(Np_values)
    Np = Np_values(i);
    win = window(window_type, Np);
    
    if Np < N
        % y_windowed = y(1:Np) .* win(1:Np);
        subplot(2, 3, i+1)
        periodo(y_windowed(1:Np), Np)
        title(['Periodogram with Np = ' num2str(Np)]);
    else
        % y_windowed = [y; zeros(Np - N, 1)] .* win;
        subplot(2, 3, i+1)
        periodo(y_windowed, Np)
        title(['Periodogram with Np = ' num2str(Np)]);
    end
end
sgtitle('Effect of different values of N on the periodogram');


%--------------------------- point C ------------------------------------
N = length(y);
y_reverse = y(N:-1:1);
convolution = conv(y, y_reverse)/numel(y);
correlation_sequence = zeros(1, N);
for i=1:N
    correlation_sequence(i) = convolution(N+i-1);
end
% Display the correlation sequence
figure;
stem(0:N-1, correlation_sequence);
title('Correlation Sequence');
xlabel('n');
ylabel('Correlation Value');
grid on;

%--------------------------- point D ------------------------------------
[pred_error_var, pred_error_filter, parcors_seq] = my_levinson(correlation_sequence);
N = 1024;
filter_f = fft(pred_error_filter, N);
filter_f = filter_f(1:(N/2)+1);
f = 0:1/N:0.5;
AR_spectrum = pred_error_var(end)./(abs(filter_f).^2);

figure;

subplot(4, 1, 1);
plot(f, 10*log10(AR_spectrum))
hold on;
periodo(y, 1024)
title('Periodogram and Autoregressive Spectrum Estimate')
legend('Autoregressive Spectrum Estimate', 'Periodogram')
xlabel('Normalized Frequency')
ylabel('dB')

subplot(4, 1, 2);
stem(0:20, pred_error_filter)
axis([0, 20, -inf, inf])
title('Filter Coefficients')
xlabel('Coefficient Index')
ylabel('A(0:20)')
grid on

subplot(4, 1, 3);
stem(0:20, pred_error_var)
axis([0, 20, -inf, inf])
title('Evolution of the Prediction Error Variance')
xlabel('Order')
ylabel('sigma^2(0:20)')
grid on

subplot(4, 1, 4);
stem(1:20, parcors_seq)
axis([1, 20, -inf, inf])
title('PARCORS')
xlabel('Order')
ylabel('K(1:20)')
grid on