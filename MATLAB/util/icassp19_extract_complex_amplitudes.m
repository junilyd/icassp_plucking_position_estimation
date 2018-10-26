function [amplitudesAbs, f0_iHS, inharmonicityCoeff, estimatedFret, estimatedString, amplitudes ] = icassp19_extract_complex_amplitudes(x, EStringTuning, fs, BModel)
  	if nargin < 4
  		BModel = 2e-4; %random number
  	end
    % if no training data uncomment these search areas.
    %betaCoeffSearchArea = [1e-5:1e-5:8e-4];
    f0Area = [75:675];
    nFFT = 2^19;
    L_HS = 5;
    betaRes = 1e-7;
    
    x = smc_apply_gaussian_window(x);
    [fAxis,X,XplotdB] = icassp19_fft(x, fs, nFFT);
    f0_HS = icassp19_harmonic_summation(X, f0Area, L_HS, fs);

    % Sorting by f0-candidates, using training data from file "mats/betaMean_Firebird.mat"
    [stringCandidates, fretCandidates, numCandidates] = smc_find_f0_candidates(f0_HS, EStringTuning);
       
    f0Model = EStringTuning*[1 2^(5/12) 2^(10/12) 2^(15/12) 2^(19/12) 2^(24/12)];
    f0Model = f0Model.*(2.^((0:12)/12))'

    mu = [f0Model(:) BModel(:)];
    muNorm = [mu(:,1)./max(abs(mu(:,1))) mu(:,2)./max(abs(mu(:,2)))]

    for nc = 1:1%numCandidates
        clear betaCoeffSearchArea;       
        % betaSearchArea by knowing string and fret
        %center = BModel(fretCandidates+1, stringCandidates);
        betaCoeffSearchArea = [1e-5:betaRes:6e-4];%[min(min(BModel)):1e-7:max(max(BModel))];% 
        %betaCoeffSearchArea =[center*(1-0.02):1e-7:center*(1+0.02)]; 

        %[X] = smc_equalize_harmonics(X, fs, f0_HS);

		L_iHS = 26;%23+floor(2500/f0_HS);
		L_classify = [1:26];
        [est_f0, betaCoeff] = icassp19_inharmonic_summation(X, f0_HS, L_iHS, fs, betaCoeffSearchArea)
    	[indeces] = icassp19_inharmonic_index(X, fs, L_iHS, est_f0, betaCoeff);
    	%costFunction   =  abs(betaCoeff-BModel(fretCandidates+1, stringCandidates)).^2  ;%sum(X(indeces(L_classify)));
    
        xMu = [est_f0./max(abs(mu(:,1))) betaCoeff./max(abs(mu(:,2)))]

        costFunction   =  sqrt( (muNorm(:,1)-xMu(:,1)).^2 + (muNorm(:,2)-xMu(:,2)).^2 );%sum(X(indeces(L_classify)));

    end
    %% Classify string and fret as ML of produced inharmonicity cost function
    [C,I] = min(costFunction);% max(costFunction);

    estimatedFret = mod(I,13)-1
    estimatedString = ceil((I+13)/13)-1
    
% 	estimatedString = stringCandidates(I);
% 	estimatedFret = fretCandidates(I);
f0_iHS = est_f0;
inharmonicityCoeff = betaCoeff;

    Z = smc_Z_inharmonic(f0_iHS,length(x),fs,L_iHS,inharmonicityCoeff); % L_iHs was 13.
    amplitudes = inv(Z'*Z)*Z'*x;
    amplitudesAbs = abs(amplitudes);
end
