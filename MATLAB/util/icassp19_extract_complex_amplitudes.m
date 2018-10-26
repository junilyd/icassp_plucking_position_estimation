function [amplitudesAbs, pitchEstimate, BEstimate, estimatedFret, estimatedString ] = icassp19_extract_complex_amplitudes(x, EStringTuning, fs, BModel)
  	if nargin < 4
  		BModel = 2e-4; %random number
  	end
    % if no training data uncomment these search areas.
    %betaCoeffSearchArea = [1e-5:1e-5:8e-4];
    f0Limits = [75 700];
    nFFT = 2^19;
    L_HS = 5;
    betaRes = 1e-7;
    
    x = smc_apply_gaussian_window(x);
    [~,X] = icassp19_fft(x, fs, nFFT);
    f0Harmonic = icassp19_harmonic_summation(X, f0Limits, L_HS, fs);
       
    f0Model = EStringTuning*[1 2^(5/12) 2^(10/12) 2^(15/12) 2^(19/12) 2^(24/12)];
    f0Model = f0Model.*(2.^((0:12)/12))';
    mu = [f0Model(:) BModel(:)];
    muNorm = [mu(:,1)./max(abs(mu(:,1))) mu(:,2)./max(abs(mu(:,2)))];

    BSearchGrid = [1e-5:betaRes:6e-4];%[min(min(BModel)):1e-7:max(max(BModel))];% 

    M = 26; %23+floor(2500/f0_HS);
    [pitchEstimate, BEstimate] = icassp19_inharmonic_summation(X, f0Harmonic, M, fs, BSearchGrid);
    [indeces] = icassp19_inharmonic_index(X, fs, M, pitchEstimate, BEstimate);
    
    xMu = [pitchEstimate./max(abs(mu(:,1))) BEstimate./max(abs(mu(:,2)))];
    euclideanDistance   =  sqrt( (muNorm(:,1)-xMu(:,1)).^2 + (muNorm(:,2)-xMu(:,2)).^2 );

    %% Classify string and fret as ML of produced inharmonicity cost function
    [C,I] = min(euclideanDistance);% max(costFunction);

    estimatedFret = mod(I,13)-1
    estimatedString = ceil((I+13)/13)-1
    
    pitchEstimate = pitchEstimate;
    BEstimate = BEstimate;

    Z = smc_Z_inharmonic(pitchEstimate,length(x),fs,M,BEstimate);
    amplitudes = inv(Z'*Z)*Z'*x;
    amplitudesAbs = abs(amplitudes);
end
