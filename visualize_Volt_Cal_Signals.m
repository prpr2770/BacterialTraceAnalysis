% Script to visualize the Volt-Calcium signals.


% =========================================================================
% 1. Read the original data:
close all; clear all;
loadData=0;% comment below line after  1 run.

if loadData == 0
    try
        inDIR = 'H:\KraljLab\studyDenoising\';
        fname = strcat(inDIR,'chg_det_volt_ca_sigs.mat');
        load(fname); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

outDIR = 'H:\KraljLab\studyDenoising\chgScores\';
% =========================================================================
% 2. Iterate over each trace and plot the volt-ca-signals

for idx = 1:length(ca_sigs_noisy)
    
    % Voltage signals
    volt_sig = volt_sigs_noisy(idx).sig;
    volt_sig_dn = volt_sigs_noisy(idx).sig_dn_wav;
    volt_sig_chg = volt_sigs_noisy(idx).sig_chg_score;
    volt_sig_dn_chg = volt_sigs_noisy(idx).sig_dn_chg_score ;
    
    % Calcium signals
    ca_sig = ca_sigs_noisy(idx).sig;
    ca_sig_dn = ca_sigs_noisy(idx).sig_dn_wav;
    ca_sig_chg = ca_sigs_noisy(idx).sig_chg_score;
    ca_sig_dn_chg = ca_sigs_noisy(idx).sig_dn_chg_score ;
    
    
    % -------------------------------------------------------------------------
    % Fig1 = Subplot(2,2)
    % 1. Volt : noisy -vs- denoisy
    % 2. Volt events  : noisy-vs-dn
    % 3. Ca : noisy-vs-dn
    % 4. Ca events : noisy-vs-dn
    
    fig1 = figure(1)
    subplot(2,2,1)
    p1 = plot(volt_sig,'g'); hold on;  
    p2 = plot(volt_sig_dn,'r'); hold off;
    legend([p1,p2],'nsy','de-nsy')
    ylabel('volt')
    subplot(2,2,2)
    plot(volt_sig_chg/max(volt_sig_chg),'g'); hold on;  plot(volt_sig_dn_chg/max(volt_sig_dn_chg),'r'); hold off;
    ylabel('chg-score')
    subplot(2,2,3)
    plot(ca_sig,'g'); hold on;  plot(ca_sig_dn,'r'); hold off;
    ylabel('calcium')
    subplot(2,2,4)
    plot(ca_sig_chg/max(ca_sig_chg),'g'); hold on;  plot(ca_sig_dn_chg/max(ca_sig_dn_chg),'r'); hold off;
    ylabel('chg-score')
    
    set(gcf,'NextPlot','add');
axes;
heading = sprintf('Volt-Ca signals for trackID:%d',idx);
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');
    
plt_nm = sprintf('volt_ca_trackID_%d.png',idx);
        saveas(fig1, fullfile(outDIR, plt_nm), 'png');
        
    % -------------------------------------------------------------------------
    % Fig2 = Subplot(2,1)
    % 1. volt_dn -vs- ca_dn
    % 2. evt_volt_dn -vs- evt_ca_dn
    
    fig2 = figure(2)
    subplot(2,1,1)
    v1 = plot(volt_sig_dn,'g'); hold on; 
    c1 = plot(ca_sig_dn,'r'); hold off;
    legend([v1,c1],'volt','calc')
    ylabel('dn-signal')
    subplot(2,1,2)
    plot(volt_sig_dn_chg/max(volt_sig_dn_chg),'g'); hold on; plot(ca_sig_dn_chg/max(ca_sig_dn_chg),'r'); hold off;
    ylabel('chg-score')
        set(gcf,'NextPlot','add');
axes;
heading = sprintf('Volt-Ca signals for trackID:%d',idx);
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');

plt_nm = sprintf('compare_volt_ca_trackID_%d.png',idx);
        saveas(fig2, fullfile(outDIR, plt_nm), 'png');
        
    % =========================================================================
    % 3. Iterate over each trace and plot the volt-ca-signals
    
    close all;
end