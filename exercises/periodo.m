  function periodo(signal,N)
% function periodo(signal,N)
% Plots the periodogram of the signal over N frequential bins.
% The spectrum (in dB) is represented over [0 0.5]*fs where fs=1  
% is the normalized sampling rate.
% N must be even.

ns=length(signal);
if rem(N,2)~=0
  error('N must be even')
end

signal=signal(:)';	% force a row vector
if N<=ns
  in=signal(1:N);
  ns=N;
else
  in=[signal zeros(1,N-ns)];	% zero padding
end
dsp=abs(fft(in));
DSP=dsp(1:(N/2)+1)/sqrt(ns);

f=0:1/N:0.5;
plot(f,10*log10(DSP.^2))
grid
xlabel('Normalized frequency')
ylabel('dB')



