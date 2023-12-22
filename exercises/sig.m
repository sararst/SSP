function [y,sigpar]=sig(N,sigpar)
%function [y,sigpar]=sig(N,sigpar)
%
%	generates N samples of a signal y, consisting of
%	two sinusoids (with amplitudes A1, A2, frequencies f1, f2,
%	and random phases) plus white noise of unit variance.
%	sigpar = vector of 4 signal parameters (f1,f2,A1,A2)
%	if [y,sigpar]=sig(N) is used then
%	sigpar=(0.057,0.082,20,20)

if nargin==1
  sigpar(1)=0.057;
  sigpar(2)=0.082;
  sigpar(3)=20;
  sigpar(4)=sigpar(3);
end;
f1=sigpar(1);
f2=sigpar(2);
A1=sigpar(3);
A2=sigpar(4);
y=A1*cos(2*pi*(f1*[1:N]'+rand))+A2*cos(2*pi*((f2)*[1:N]'+rand));
y=y+randn(N,1);
y=y(:);
