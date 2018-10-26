% -------------------------------------------------------------------
% calcualates the FFT of the input signal. returns up to nyquist freq.
% It does the zeropadding up to next power of 2 
% with the built a MATLAB function.
%
% INPUT: 
%       sig:  signal to analyse
%       Fs:   sample rate of that signal
%       N     length of the FFT, which can change the resolution of the plot
% OUTPUT:
%       f_axis:     x-axis for plotting
%       fft_sig:    linear fft signal
%       fft_sig_dB: logarithmic fft signal in dB
%       NFFT:       the power of 2 used for the fft
%      
% EXAMPLE: 
%         [sig,fs] = audioread('testfiles/A_mid.wav');
%         [f,fft_sig, fft_dB] = smc_fft(sig,fs);
%         plot(f,fft_sig); figure;
%         plot(f,fft_dB);
%         
% -----------------------------------------------------------
% function [f_axis, fft_sig, fft_sig_dB] = fft_(sig, Fs, N)    
% -----------------------------------------------------------
function [f_axis, fft_sig, fft_sig_dB] = icassp19_fft(sig,Fs,N)    
    fft_long = fft(sig,N);
    fft_sig  = 2*abs(fft_long(1:N/2)).^2;
    fft_sig_dB = 10*log10(fft_sig);
    f_axis   = Fs/2 * linspace(0,1,N/2);
end
