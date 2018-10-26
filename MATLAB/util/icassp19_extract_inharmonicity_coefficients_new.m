% -------------------------------------------------------------
% Extracts the inharmonicity coefficient (betacoeff)
% it also estimates string and fret position for a 
% standard guitar tuning. 
% 
%   INPUTS:     
%               x              : Hilbert transform of audio signal.
%               EStringTuning  : pitch estimate of lowest string.
%                                (assuming standard guitar tuning.)
%               fs             : sample rate of the input [x].
%               betaMean       : statistical reference table (optional). 
%               L              : The desired number of harmonics
%
%   OUTPUTS:
%               estimatedString: String estimate based on EStringTuning
%               estimatedFret  : Fret estimate based on EStringTuning
%               betaCoeff      : Inharmonicity coefficient
%               f0_HS          : course pitch estimate
%               f0_iHS         : inharmonic pitch estimate
%               fAxis          : f axis for plotting th fft
%               XplotdB        : real fft vector
%               diffMedianMean : difference between binMean and binMedian vectors.
%               binMeanVector  : mean vector, containg mean value across each window
%               XOriginal      : original fft vector as output
%               binMedianVector: median vector, containg mean value across each window
%               threshVector   : threshhold vector at 90% above mean level, relative to peak level. 
% -------------------------------------------------------------------------
% [estimatedString, estimatedFret, betaCoeff, f0_HS, f0_iHS, fAxis, XplotdB, ...
% diffMedianMean, binMeanVector, XOriginal, binMedianVector, ThreshVector ] ... 
% = smc_extract_inharmonicity_coefficients(x, EStringTuning, fs, betaMean)
% ------------------------------------------------------------------------
%
%

%
function [estimatedString, estimatedFret, betaCoeff, f0_HS, f0_iHS, fAxis, XplotdB, X ] = icassp19_extract_inharmonicity_coefficients_new(x, EStringTuning, fs, BMean)
  	if nargin < 4
  		BMean = 2e-4; %random number
  	end
    % if no training data uncomment these search areas.
    %betaCoeffSearchArea = [1e-5:1e-5:8e-4];
    f0Area = [75:675];
    nFFT = 2^19;
    L_HS = 5;
    
    %x = smc_apply_gaussian_window(x);
    [fAxis,X,XplotdB] = icassp19_fft(x, fs, nFFT);
    f0_HS = icassp19_harmonic_summation(X, f0Area, L_HS, fs);

    % Sorting by f0-candidates, using training data from file "mats/betaMean_Firebird.mat"
    [stringCandidates, fretCandidates, numCandidates] = smc_find_f0_candidates(f0_HS, EStringTuning);

    for nc = 1:numCandidates
        clear betaCoeffSearchArea;       
        % betaSearchArea by knowing string and fret
        center = BMean(fretCandidates(nc)+1, stringCandidates(nc));
        X=X;
		L_iHS = 23+floor(2500/f0_HS);
		L_classify = [1:26];
        [est_f0(nc), betaCoeff(nc)] = smc_inharmonic_summation_tuner(X, f0_HS, L_iHS, fs, betaCoeffSearchArea);
    	[indeces] = icassp19_inharmonic_index(X, fs, L_iHS, est_f0(nc), betaCoeff(nc));
    	costFunction(nc)   = sum(X(indeces(L_classify)));
    end
    %% Classify string and fret as ML of produced inharmonicity cost function
    [C,I] = max(costFunction);
	estimatedString = stringCandidates(I);
	estimatedFret = fretCandidates(I);
	f0_iHS = est_f0(I);
end