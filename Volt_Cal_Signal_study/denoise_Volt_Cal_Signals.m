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

cellFamily = 'styph'; % 'styph'

if loadData == 0
    try
        tracesHome = 'H:\KraljLab\studyDenoising\twoSpecies_transferFunction\SecondRun';
        tracesDirName = [tracesHome filesep cellFamily];
        

        %{
        switch cellFamily
            case 'ecoli'
                tracesFileName = strcat(tracesDirName,'ecoli_data_voltage_calcium.mat');
            case 'styph'
                tracesFileName = strcat(tracesDirName,'styph_data_voltage_calcium.mat');
        end
        %}
        
        %         tracesFileName = strcat(tracesDirName,cellFamily,'volt_ca_data.mat');
        tracesFileName = [tracesDirName filesep 'volt_ca_data.mat'];
        
        load(tracesFileName); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

% outDIR = [tracesDirName filesep cellFamily];
outDIR = tracesDirName;
% =========================================================================
% 2. Extract data into structures and save.
% volt_sigs_noisy; % empty cell structures
% ca_sigs_noisy;

totalTraces = length(data);

file_nm = [cellFamily '_' 'volt_ca_sigs.mat'];
matfile_nm = fullfile(outDIR, file_nm);
matfile_obj = matfile(matfile_nm, 'Writable', true);

% intialize the data-structures
ca_sigs_noisy(1).sig = [];
ca_sigs_noisy(1).sig_dn_wav = [];
ca_sigs_noisy(1).err_dn_wav = [];
volt_sigs_noisy(1).sig = [];
volt_sigs_noisy(1).sig_dn_wav = [];
volt_sigs_noisy(1).err_dn_wav = [];



for idx = 1:totalTraces
    volt = data(idx).volt;
    ca = data(idx).calcium;
    
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

save(matfile_nm,'ca_sigs_noisy', 'volt_sigs_noisy','-v7.3');



%{
% =========================================================================
% Event Change Detection on the Noisy and Denoised Signal


for idx = 1:totalTraces
    sprintf('trace %d out of %d \n',idx,totalTraces)
    % Template for computing/updating change-detection
    %     sig = volt_sigs_noisy(idx).sig;
    %     chg_score = getChangeDetectionScore(sig);
    %     volt_sigs_noisy(idx).sig_chg_score = chg_score;
    
    
    % Voltage Signals:
    chg_volt_sigs_noisy(idx).sig_chg_score = getChangeDetectionScore(volt_sigs_noisy(idx).sig);
    chg_volt_sigs_noisy(idx).sig_dn_chg_score = getChangeDetectionScore(volt_sigs_noisy(idx).sig_dn_wav);
    
    % Calcium Signals:
    chg_ca_sigs_noisy(idx).sig_chg_score = getChangeDetectionScore(ca_sigs_noisy(idx).sig);
    chg_ca_sigs_noisy(idx).sig_dn_chg_score = getChangeDetectionScore(ca_sigs_noisy(idx).sig_dn_wav);
    
    
    % --------------------------------------------------
    % Extract peak-interval distribution for the denoised signals.
    % obtain the maxValue and time instant
    minPeakHeight = 0.6;
    minPeakProminence = 0.1;
    leftMargin = 60;
    
    % Voltage Denoised Signals:
    sig = chg_volt_sigs_noisy(idx).sig_dn_chg_score;
    [ pks, locs] = getPeakLocs_from_ChangeScore(sig,leftMargin,minPeakHeight, minPeakProminence);
    chg_volt_sigs_noisy(idx).sig_dn_chg_score_pk_locs = locs;
    
    % Calcium Denoised Signals:
    sig = chg_ca_sigs_noisy(idx).sig_dn_chg_score;
    [ pks, locs] = getPeakLocs_from_ChangeScore(sig,leftMargin,minPeakHeight, minPeakProminence);
    chg_ca_sigs_noisy(idx).sig_dn_chg_score_pk_locs = locs;
    
end

file_nm = [cellFamily '_' 'chg_det_volt_ca_sigs.mat'];
save(fullfile(outDIR, file_nm),'chg_ca_sigs_noisy', 'chg_volt_sigs_noisy');

%}