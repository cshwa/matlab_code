% FFTPLOT(in)
% This function plots the power spectrum on the input vector.

function out = fftplot(in)

a = fft(in);
a(1) = [];
b = length(a)/2;
power = abs(a(1:b)).^2
nyq = 1/2;
freq = (1:b)/b*nyq
plot(freq,power);
out = 'Done';