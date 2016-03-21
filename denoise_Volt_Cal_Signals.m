% % 12 March: Touri's tasks:
% 1. Denoising the traces using two methods: Wavelets
% 2. Determine error between the noisy and denoised signals.
% 3. Determine statistics of the error signals
% 4. Extract change-score for the denoised signals
% 5. Extract significant event-chg-locations. 


% =========================================================================
% 1. Read the original data:
close all; clear all;

loadData=0;% comment below line after  1 run.

if loadData == 0
    try
        tracesDirName = 'H:\KraljLab\';
        tracesFileName = strcat(tracesDirName,'data_voltage_calcium.mat');
        load(tracesFileName); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

outDIR = 'H:\KraljLab\studyDenoising\';
% =========================================================================
% 2. Extract data into structures and save.
% volt_sigs_noisy; % empty cell structures
% ca_sigs_noisy;

for idx = 1:length(data)
    volt = data(idx).TimeTraceVolt;
    ca = data(idx).TimeTraceCa;
    
    % remove DC component from signal
    volt  = volt - mean(volt);
    ca = ca - mean(ca);
    % archive it
    volt_sigs_noisy(idx).sig = volt ;
    ca_sigs_noisy(idx).sig = ca ;
    
    % extract wavelet denoising
    volt_denoise_wav = getWaveletDenoisedTrace(volt);
    ca_denoise_wav = getWaveletDenoisedTrace(ca);
    %archive it.
    volt_sigs_noisy(idx).sig_dn_wav = volt_denoise_wav  ;
    ca_sigs_noisy(idx).sig_dn_wav = ca_denoise_wav  ;
    
    % determine the error obtained from the wavelet denoising
    volt_err = volt_denoise_wav - volt;
    ca_err = ca_denoise_wav - ca;
    volt_sigs_noisy(idx).err_dn_wav = volt_err  ;
    ca_sigs_noisy(idx).err_dn_wav = ca_err  ;
    
    
end

file_nm = 'volt_ca_sigs.mat';
save(fullfile(outDIR, file_nm),'ca_sigs_noisy', 'volt_sigs_noisy');

% =========================================================================
% Event Change Detection on the Noisy and Denoised Signal


for idx = 1:length(data)
    % Template for computing/updating change-detection
    %     sig = volt_sigs_noisy(idx).sig;
    %     chg_score = getChangeDetectionScore(sig);
    %     volt_sigs_noisy(idx).sig_chg_score = chg_score;
    
    
    % Voltage Signals:
    volt_sigs_noisy(idx).sig_chg_score = getChangeDetectionScore(volt_sigs_noisy(idx).sig);
    volt_sigs_noisy(idx).sig_dn_chg_score = getChangeDetectionScore(volt_sigs_noisy(idx).sig_dn_wav);
    
    % Calcium Signals:
    ca_sigs_noisy(idx).sig_chg_score = getChangeDetectionScore(ca_sigs_noisy(idx).sig);
    ca_sigs_noisy(idx).sig_dn_chg_score = getChangeDetectionScore(ca_sigs_noisy(idx).sig_dn_wav);
    
    
    % --------------------------------------------------
    % Extract peak-interval distribution for the denoised signals.
    % obtain the maxValue and time instant
    minPeakHeight = 0.6;
    minPeakProminence = 0.1;
    leftMargin = 60;
    
    % Voltage Denoised Signals:
    sig = volt_sigs_noisy(idx).sig_dn_chg_score;
    [ pks, locs] = getPeakLocs_from_ChangeScore(sig,leftMargin,minPeakHeight, minPeakProminence)
    volt_sigs_noisy(idx).sig_dn_chg_score_pk_locs = locs;
    
    % Calcium Denoised Signals:
    sig = ca_sigs_noisy(idx).sig_dn_chg_score;
    [ pks, locs] = getPeakLocs_from_ChangeScore(sig,leftMargin,minPeakHeight, minPeakProminence)
    ca_sigs_noisy(idx).sig_dn_chg_score_pk_locs = locs;
    
end

file_nm = 'chg_det_volt_ca_sigs.mat';
save(fullfile(outDIR, file_nm),'ca_sigs_noisy', 'volt_sigs_noisy');