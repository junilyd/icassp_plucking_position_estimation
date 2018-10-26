% cd ~/Documents/MATLAB/code/string_classification/inharmonicity/
% addpath ~/repositories/guitar_string_finger_and_pluck_estimation/
% 
clear all;

segmentDuration = 40e-3; % ms

addpath util

%load ~/repositories/guitar_string_finger_and_pluck_estimation/mats/betaMeanfirebird_40ms_from_12th_fret_betares0.10u_nFFT2^19.mat
%load trained_model_from_12th_fret.mat;

load B12.mat;
% 
%f0Area = [65:0.1:700];
%L_HS = 5;
mirverbose(0);

[recordedSignal,fs]=audioread('~/repositories/guitar_string_finger_and_pluck_estimation/model_based_estimator_paper/audio_recordings/tablature_constant_note22.wav');

EStringTuning = 83.17;%mu(1,1,2).*normalizationFactor(2);

[segment, onsetInSeconds] = icassp19_segment_from_all_onsets(recordedSignal,fs,segmentDuration); 
% for n = 1:size(segment,2)
% 	clear x; clear X; clear f_axis;
% 	[x] = icassp19_hilbert_transform(segment(:,n),fs); 
% [string(n), fret(n), betaCoeff, f0_HS, f0_iHS, fAxis, XplotdB, XOriginal ] ...
%     = icassp19_extract_inharmonicity_coefficients(x, EStringTuning, fs, mu(:,:,1).*normalizationFactor(1)); 
% end

for n = 1:size(segment,2)
	clear x; clear X; clear f_axis Z amplitudes amplCompl;
    [x] = icassp19_hilbert_transform(segment(:,n),fs); 
    
    [amplitudes, f0, B, fret(n), string(n), amplCompl] ...
        = icassp19_extract_complex_amplitudes(x, EStringTuning,fs, betaModelApproximation);%mu(:,:,1).*normalizationFactor(1)); 
    amplitudes = amplitudes';


L=64.3;
h = 1;
vibratingStringLength = L * 2^(-fret(n)/12);
[pluckCmFromBridge(n), Cn] = icasssp19_plucking_position_estimator_LSD(amplitudes,vibratingStringLength);% icasssp19_plucking_position_estimator_LSD
end

%%
recordDuration = 12;
stringVector = [6:-1:6]';
timeAxis = [0:1/fs:recordDuration-1/fs];
%close all
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

text(onsetInSeconds(nn)-0.2,string(nn),sprintf('%1.0f',fret(nn)),'fontsize', 25)
        end
set(gca,'ytick',[1 2 3 4 5 6])
yticklabels({'1','2','3','4','5','6'})
ylabel('String Est.'); xlabel('Time [sec]');


