% ----------------------------------------------------------------------
% Estimation of pitch and string stiffness from using inharmonic summation
%
%   INPUTS:
%           X:            fft of input signal
%           pitchInput :  first estimate of pitch (usually from harmonic_summation_tuner())
%           L:            Number of Harmonics to use.
%           fs:           Sample rate to use
%           betaArea :    Searsh area for stiffnessCoeffs 
%   OUTPUTS:
%           pitch:        The pitch estimate
%   BstiffnessEstimate :  The stiffness coefficients estimate
% ------------------------------------------------------------------------------------------------
% [inharmonicPitch, stiffnessCoeff] = smc_inharmonic_summation_tunerv3(X, f0_area, L, fs,betaArea)
% ------------------------------------------------------------------------------------------------
function [pitchEstimate, BEstimate, costFunction, betaArea, pitchArea] = icassp19_inharmonic_summation(X, pitchInput, M, fs, betaArea, nFFT)
%pitchInput=pitch(1);
if ~exist('nFFT'), nFFT = 2^19; end

pitchInput = [pitchInput 2*pitchInput];
pitchWidth = 3*nFFT/2^18*fs/nFFT;% 3 is for the gaussian window%pitchInput(1)-pitchInput(1)*0.995;

for pp = 1:1 % testing lower(1) and upper(2) octave
    B=1;
    for beta=betaArea
        pitchArea = [pitchInput(pp)-pitchWidth:fs/nFFT:pitchInput(pp)+pitchWidth];
        for pII=1:length(pitchArea)
            [HL] = icassp19_inharmonic_index(X, fs, M, pitchArea(pII),beta);
            HL(HL>=length(X))=[];
            costFunction(:,B,pII) = sum(X(HL));% X(HL)'*((1:length(HL))'.^1);
        end
        B=B+1;
    end
    
    betaAreaSize  = size(costFunction,2);
    pitchAreaSize = size(costFunction,3);
    %% CostFunction
    [C(pp) I] = max(costFunction(:));
    
    %% Pitch
    pitchIndex = ceil(I/betaAreaSize);
    inharmonicPitchV2(pp) = pitchArea(pitchIndex);

    %% Beta
    BstiffnessIndex = mod(I,betaAreaSize);
    if BstiffnessIndex == 0, 
        BstiffnessIndex = betaAreaSize; 
    end
    BCandidates(pp) = betaArea(BstiffnessIndex);
end
%% loop end
[octaveMax octaveMaxIndex] = max(C);
pitchEstimate = inharmonicPitchV2(octaveMaxIndex);
BEstimate = BCandidates(octaveMaxIndex);
end
