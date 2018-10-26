% ----------------------------------------------------------------------
% Calculate analytic complex signal from real signal   
%   Calcualates the analytic signal from a real signal into size of next
%   power of 2 from lenght(sig).
%   
%   First the function calculates the zero padded FFT. Next step is to make
%   a vector, h, with the length of FFT--> NFFT. On first place in the vector 
%   a 1 is loaded. Also in NFFT/2+1 place a 1 is loaded. 
%   From h(2):NFFT/2 a 2 is loaded. 
%   A new vector is made, with same size as the previous h. The vector 
%   is loaded with the FFT result multiplied the h vector. This means the
%   last half will be zero and the first part is doubled in magnitude,
%   which result in same signal magnitude. 
%   
%
%   INPUT: 
%           sig: real signal
%   OUTPUT:
%           x: analytic signal (complex signal)
% 	    X: FFT signal
%      f_axis: values for f axis (mostly for ploting cases)
%
% --------------------
% x = smc_hilbert(sig)
% --------------------
function [x, X, f_axis] = icassp19_hilbert_transform(sig, fs)
if (~exist('fs')), fs = 44.1e3; end;

N = length(sig);
if mod(length(sig),2) == 1
   sig = sig(1:end-1);
   N = N-1;
end
fft_long = fft(sig, N);
X = 2*abs(fft_long(1:N/2)).^2; 
f_axis   = fs/2 * linspace(0,1,N/2);
h = zeros(1,N)';
h(1) = 1;
h(2:N/2) = 2;
h(N/2+1) = 1;
x_fft = zeros(1,N)';
for i = 1:N
    x_fft(i) = fft_long(i)*h(i);
end
x = ifft(x_fft);
end
