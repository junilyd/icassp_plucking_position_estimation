% cd ~/Documents/MATLAB/code/string_classification/inhar-3monicity/
% addpath ~/repositories/guitar_string_finger_and_pluck_estimation/
% 
clear all;
mirverbose(0);
addpath util

%% Load trained model (training data was captured from the 12th fret)

%load ~/repositories/guitar_string_finger_and_pluck_estimation/mats/betaMeanfirebird_40ms_from_12th_fret_betares0.10u_nFFT2^19.mat
%load trained_model_from_12th_fret.mat;
load B12.mat;
BModel = betaModelApproximation;
EStringTuning = 83.17;%mu(1,1,2).*normalizationFactor(2);
f0Model = EStringTuning*[1 2^(5/12) 2^(10/12) 2^(15/12) 2^(19/12) 2^(24/12)];
f0Model = f0Model.*(2.^((0:12)/12))';
mu = [f0Model(:) BModel(:)];
muNorm = [mu(:,1)./max(abs(mu(:,1))) mu(:,2)./max(abs(mu(:,2)))];

%% Initialize implementation constants 
segmentDuration = 40e-3; % segment duration in seconds.
LOpen = 64.3; % assumed length of all open strings.
M = 25; % assumed high number of harmonics (M>>1).
MInitial = 5; % number of harmonics for initial harmonic estimate (B=0).
f0Limits = [75 700]; % boundaries for f0 search grid in Hz.
nFFT = 2^19; % Length of  zero-padded FFT.
betaRes = 1e-7; % resolution of search grid for B in m^(-2)
BSearchGrid = [min(min(BModel)):1e-7:max(max(BModel))]; % searh grid for B. 


%% read in the observed guitar recording and do onset detection
[recordedSignal,fs]=audioread('~/repositories/guitar_string_finger_and_pluck_estimation/model_based_estimator_paper/audio_recordings/tablature_constant_note22.wav');
% segment the signal from every onset event (i.e. 40 ms)
[segment, onsetInSeconds] = icassp19_segment_from_all_onsets(recordedSignal,fs,segmentDuration); 

%% Estimate string, fret and plucking position on every 40 ms. segment
for n = 1:size(segment,2)
    % Hilbert transform and windowing
    [x] = icassp19_hilbert_transform(segment(:,n),fs); 
    x = icassp19_apply_gaussian_window(x);

    %% Feature extraction with the inharmonic pitch estimator
    % The implementation of Eq. (17) is done with one FFT, since it is fast
    % Hence, it is equivalent to harmonic summation which we extend to
    % inharmonic summation.
    % See details on harmonic summation in Christensen and Jakobsson [27].
    [~,X] = icassp19_fft(x, fs, nFFT);
    f0Harmonic = icassp19_harmonic_summation(X, f0Limits, MInitial, fs);
    [pitchEstimate, BEstimate] = icassp19_inharmonic_summation(X, f0Harmonic, M, fs, BSearchGrid);
    
    % feature vector computed from the observation and normalized for
    % euclidean distance. We use the trained model as part of normalization.
    w = [pitchEstimate./max(abs(mu(:,1))) BEstimate./max(abs(mu(:,2)))];
  
    %% Classifation of String and Fret (maximum likelihood w. uniform prior)
    euclideanDistance   =  sqrt( (w(:,1)-muNorm(:,1)).^2 + (w(:,2)-muNorm(:,2)).^2 );
    [C,I] = min(euclideanDistance);
    fretEstimate(n) = mod(I,13)-1; % <-- due to a matrix structure (13x6)
    stringEstimate(n) = ceil((I+13)/13)-1;

    %% Estimate the amplitudes (alpha vector)
    Z = smc_Z_inharmonic(pitchEstimate,length(x),fs,M,BEstimate);
    alpha = inv(Z'*Z)*Z'*x;
    amplitudesAbs = abs(alpha)'; % absolute values for the estaimtor

    %% Plucking Position Estimator (minimizer of log spectral distance)
    L = LOpen * 2^(-fretEstimate(n)/12);
    [pluckCmFromBridge(n)] = icasssp19_plucking_position_estimator_LSD(amplitudesAbs,L);% icasssp19_plucking_position_estimator_LSD
end

%% Plot the results
recordDuration = 12;
timeAxis = [0:1/fs:recordDuration-1/fs];
figure(10);
    subaxis(3,1,1,'sh',.1,'sv',.3,'mt',0,'pt',0.01);

        plot(timeAxis,recordedSignal); ylabel('Ampl.');
        set(gca,'xticklabel',[]);
    subaxis(3,1,1,'sh',.3,'sv',.01,'pt',0.03,'margin',0.1);

        plot(onsetInSeconds(1:n),pluckCmFromBridge(1:n),'x');
        ylim([4 32])
        set(gca,'xticklabel',[]); ylabel('$\hat{P}$[cm]','interpreter','latex')
        grid minor;
    subaxis(3,1,2,'sh',.03,'sv',0,'mt',0,'pt',0.09);

        for p=1:6
            plot([0,recordDuration],[p,p], 'Color', [0.4 0.4 0.4], 'linewidth',1); hold on;
        end
        ylim([0.1 6.9]);
        for nn = 1:n

text(onsetInSeconds(nn)-0.2,stringEstimate(nn),sprintf('%1.0f',fretEstimate(nn)),'fontsize', 20)
        end
set(gca,'ytick',[1 2 3 4 5 6])
yticklabels({'1','2','3','4','5','6'})
ylabel('String Est.'); xlabel('Time [sec]');


