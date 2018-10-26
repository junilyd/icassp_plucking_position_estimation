%   ----------------------------------------------------------
%   Estimates the pitch from harmonic summation
%   Here the tuner also zooms in to get some calc. efficiency
%   --------------------------------------------------------------------------------
%   INPUTS:
%           X:       fft of input signal
%           f0_area: Vector of fundamental frequrncies in hertz (the search area)
%           L:       Number of Harmonics to use.
%           fs:      Sample rate to use (usually samplerate of input signal)
%   OUTPUTS:
%           pitch:   The pitch estimate
% --------------------------------------------------------
% pitch3 = smc_harmonic_summation_tuner(X, f0_area, L, fs)
% --------------------------------------------------------
function [pitch3,f0Index] = icassp19_harmonic_summation(X, f0_area, L, fs)
if nargin < 5
f0_area = [min(f0_area):0.1:max(f0_area)];
end


i=1;
for f=f0_area
    [index] = icassp19_harmonic_index(X, fs, L, f);
    cost(i) = sum(X(index)); i = i+1;
end
[C, I] = max(cost);
pitch = f0_area(I);

f0_area2 = [pitch-4:0.01:pitch+4]; % changed from +-4 Hz
i=1;
for f=f0_area2
    [index] = icassp19_harmonic_index(X, fs, L, f);
    cost2(i) = sum(X(index));i=i+1;
end
[C,I] = max(cost2);
pitch2 = f0_area2(I);

f0_area3 = [pitch2-0.01:0.001:pitch2+0.01];
i=1;
for f=f0_area3
    [index] = icassp19_harmonic_index(X, fs, L, f);
    cost3(i) = sum(X(index));i=i+1;
end
[C,I] = max(cost3);
pitch3 = f0_area3(I);

% return maximum near f0 estimate
[firstPeakIndex] = icassp19_harmonic_index(X, fs, 1, pitch3);
[lowerIndex] = icassp19_harmonic_index(X, fs, 1, pitch3-pitch3/4);
[upperIndex] = icassp19_harmonic_index(X, fs, 1, pitch3+pitch3/4);

[~,peakNdx] = max(X(lowerIndex:upperIndex));
f0Index=peakNdx+lowerIndex-1;


end
