% Apply Sinc Filter to the signals
% Reference: https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
close all; clear all;


Fs = 5;
fc = 0.5; % Cutoff frequency as a fraction of the sampling rate (in (0, 0.5)).
b = 0.08; % Transition band, as a fraction of the sampling rate (in (0, 0.5)).
N = 64; % ceil(4/b)  >> How to incorporate this?
n = 0:1:N-1;
sinc_filt = sinc(2*fc*(n - (N-1)/2));
blkman_window = 0.42 - 0.5*cos(2*pi*n/(N-1))+0.08*cos(4*pi*n/(N-1));
sinc_window = blkman_window.*sinc_filt;
% sinc_window_final = 2*fc*sinc(2*fc*(n - (N-1)/2)).*blkman_window;
figure;
plot(n,sinc_window,'g',n,blkman_window,'r',n,sinc_filt,'b')

% =========================================================================
% 1. Read the original data:

loadData=0;% comment below line after  1 run.

if loadData == 0
    try
        inDIR = 'H:\KraljLab\studyDenoising\';
        fname = strcat(inDIR,'volt_ca_sigs.mat');
        load(fname); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end
% =========================================================================
close all;
for idx = 1:2%length(ca_sigs_noisy)
    
    sig = volt_sigs_noisy(idx).sig;
    clear_sig = conv(sig,sinc_window,'same');
    
    figure;
    plot(sig,'g'); hold on;
    plot(clear_sig,'r')
    
end

% =========================================================================
