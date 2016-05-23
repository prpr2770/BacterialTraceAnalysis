% Studying the transfer function relationship between the Calcium and
% Voltage signals.


% Procedure:
% 1. Read each dataset pair of volt-ca signals.
% 2. create the iddata object
% 3. obtain transfer function relationship containing one/two poles.
% (Q: How does the number of poles of the system matter? What about the number of zeros?)
% alternate the input and output relationship: volt-ca, ca-volt
% 4. What do you understand from linear systems analysis?

% -------------------------------------------------------------------------

% read all the signals.
close all;

% %{
clear all;


%%
inDIR = 'H:\KraljLab\studyDenoising\';
outDIR = [inDIR filesep 'fft_transferFunctionStudy\'];
if ~exist(outDIR)
    mkdir(outDIR)
end

caDIR = [outDIR filesep 'inputCA\'];
voltDIR = [outDIR filesep 'inputVOLT\'];

if ~exist(caDIR) || ~exist(voltDIR)
    mkdir(caDIR)
    mkdir(voltDIR)
end



%%
% Read/extract dataset.

file_nm = 'chg_det_volt_ca_sigs.mat';
load([inDIR filesep file_nm]);

% %}

inputType = 'ca'; %'volt' 'ca'


Fs = 5;
Ts = 1/Fs;



% ensure size of the variable-structures are equal.
if length(volt_sigs_noisy) ~= length(ca_sigs_noisy)
    error('Unequal dataset lengths.')
    
else
    
    for idx = 1:5%length(volt_sigs_noisy)
        volt_sig = volt_sigs_noisy(idx).sig_dn_wav';
        ca_sig = ca_sigs_noisy(idx).sig_dn_wav';
        
        
        
        switch inputType
            
            case 'ca'
                saveDIR = caDIR;
                u = ca_sig;
                y = volt_sig;
                
            case 'volt'
                saveDIR = voltDIR;
                u = volt_sig;
                y = ca_sig;
        end
        
        u_fft = fft(u);
        y_fft= fft(y);
        
        epsilon = 0.02;
        idx_nz_u_fft = find(u_fft > epsilon);
        
        
        transfer_fft = y_fft;
        transfer_fft(idx_nz_u_fft) = y_fft(idx_nz_u_fft)./u_fft(idx_nz_u_fft);
        
        abs_transfer_fft = abs(transfer_fft);
        phi_transfer_fft = angle(transfer_fft);
        
        fig1 = figure(1)
        N = length(u);
        t = Ts*(1:N);
        subplot(2,1,1)
        plot(t,u)
        subplot(2,1,2)
        plot(t,y)
        
        fig2 = figure(2)
        subplot(2,1,1)
        plot(log(abs_transfer_fft))
        
        subplot(2,1,2)
        plot(phi_transfer_fft)
        
        plt_nm = sprintf('compare_io_model_%d.png',idx);
        saveas(fig1, fullfile(saveDIR, plt_nm), 'png');
        
        pause
        
        
        %%
        close all;
    end
end