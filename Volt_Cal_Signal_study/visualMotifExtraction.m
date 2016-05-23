% Script to extract motifs from interactive plots of traces.

% SHOULD THE MOTIFS BE NORMALIZED BEFORE COMPUTING THE EIGEN-MOTIFS?

% -------------------------------------------------------------------------

% read all the signals.
close all;

% %{
clear all;


%%
inDIR = 'H:\KraljLab\studyDenoising\';
outDIR = [inDIR filesep 'visualMotifs\'];
if ~exist(outDIR)
    mkdir(outDIR)
end




%%
% Read/extract dataset.

file_nm = 'chg_det_volt_ca_sigs.mat';
load([inDIR filesep file_nm]);

% %}

Fs = 5;
Ts = 1/Fs;



% ensure size of the variable-structures are equal.
if length(volt_sigs_noisy) ~= length(ca_sigs_noisy)
    error('Unequal dataset lengths.')
    
else
    
    for idx = 1:length(volt_sigs_noisy)
        
        % volt-signals
        volt_sig = volt_sigs_noisy(idx).sig_dn_wav';
        volt_cgsc = volt_sigs_noisy(idx).sig_chg_score';
        volt_dn_cgsc = volt_sigs_noisy(idx).sig_dn_chg_score';
        
        % ca-signals
        ca_sig = ca_sigs_noisy(idx).sig_dn_wav';
        ca_cgsc = ca_sigs_noisy(idx).sig_chg_score';
        ca_dn_cgsc = ca_sigs_noisy(idx).sig_dn_chg_score';
        
        
        % plot the signal and extract motif locations.
        
        
        % voltage motifs
        
        fig1 = figure(1)
        subplot(2,1,1)
        plot(volt_sig)
        subplot(2,1,2)
        plot(volt_dn_cgsc/max(volt_dn_cgsc),'r');
        hold on;
        plot(volt_cgsc/max(volt_cgsc),'g');
        hold off;
        
        [x,y] = ginput(10)
        
        volt_motifs(idx).motif_locs = floor(x');
        
        
        % extract calcium motifs
        
        fig2 = figure(2)
        subplot(2,1,1)
        plot(ca_sig)
        subplot(2,1,2)
        plot(ca_dn_cgsc/max(ca_dn_cgsc),'r');
        hold on;
        plot(ca_cgsc/max(ca_cgsc),'g');
        hold off;
        
        [x,y] = ginput(10)
        
        ca_motifs(idx).motif_locs = floor(x');
        
        
        %%
        close all;
    end
end


% archive the motif locations
file_nm = 'visual_motifs.mat';
save(fullfile(outDIR, file_nm),'ca_motifs', 'volt_motifs');


%% Extract the motifs subsequences and then extract eigen-motifs.
all_ca_motifs = [];
all_volt_motifs = [];

ca_motif_size = 60;
volt_motif_size = 30;

for idx = 1: length(ca_motifs)
    
    volt_sig = volt_sigs_noisy(idx).sig_dn_wav;
    ca_sig = ca_sigs_noisy(idx).sig_dn_wav;
    
    % ---------------------------------------------------------
    % volt motifs
    motif_beg = volt_motifs(idx).motif_locs;
    motif_end = volt_motif_size + motif_beg;
    
    all_motifs = [];
    sig = volt_sig;
    for i = 1:length(motif_beg)
        motif = sig(motif_beg(i):motif_end(i));
        all_motifs = [all_motifs; motif];
    end
    volt_motifs(idx).motif_subseqs = all_motifs;
    all_volt_motifs = [all_volt_motifs; all_motifs];
    
    % ---------------------------------------------------------
    % ca-motifs
    motif_beg = ca_motifs(idx).motif_locs;
    motif_end = ca_motif_size + motif_beg ;
    
    all_motifs = [];
    sig = ca_sig;
    for i = 1:length(motif_beg)
        motif = sig(motif_beg(i):motif_end(i));
        all_motifs = [all_motifs; motif];
    end
    ca_motifs(idx).motif_subseqs = all_motifs;
    all_ca_motifs = [all_ca_motifs; all_motifs];
    
end

file_nm = 'visual_motifs.mat';
save(fullfile(outDIR, file_nm),'ca_motifs', 'volt_motifs');

file_nm = 'all_visual_motifs.mat';
save(fullfile(outDIR, file_nm),'all_ca_motifs', 'all_volt_motifs');

%% extract eigenMotifs

% SHOULD THE MOTIFS BE NORMALIZED BEFORE COMPUTING THE EIGEN-MOTIFS?

% voltage
source = all_volt_motifs;
target = source;
k = 12;
[S W]= GeneralizedGreedySelection(target, source, k)
volt_eigMotifs = source(S,:);

% calcium
source = all_ca_motifs;
target = source;
k = 12;
[S W]= GeneralizedGreedySelection(target, source, k)
ca_eigMotifs = source(S,:);


file_nm = 'all_visual_motifs.mat';
save(fullfile(outDIR, file_nm),'all_ca_motifs', 'all_volt_motifs','ca_eigMotifs', 'volt_eigMotifs');


% visualize eigenMotifs
%%

% visualize both simultaneously
fig5 = figure(5)
numPlots = max(size(volt_eigMotifs,1),size(ca_eigMotifs,1));
numRows = min(4, numPlots);
numCols = ceil(numPlots/numRows);

for idy =1:numPlots
    subplot(numRows,numCols,idy)
    if idy <= size(volt_eigMotifs,1)
        p1 = plot(volt_eigMotifs(idy,:),'r'); hold on;
    else
        p1 = plot(zeros(size(volt_eigMotifs(1,:))),'r'); hold on;
    end
    
    if idy <= size(ca_eigMotifs,1)
        
        p2 = plot(ca_eigMotifs(idy,:),'g'); hold off;
    else
        p2 = plot(zeros(size(ca_eigMotifs(1,:))),'r'); hold on;
    end
    if idy == 1
        legend([p1 p2],'V','Ca');
    end
    
end

plt_nm = sprintf('indiv_eigMotifs.png');
saveas(fig5, fullfile(outDIR, plt_nm), 'png');

% visualize volt-motifs
fig6 = figure(6)
subplot(2,1,1)
plot(volt_eigMotifs')
ylabel('volt')
subplot(2,1,2)
plot(ca_eigMotifs')
ylabel('ca')

plt_nm = sprintf('all_eigMotifs.png');
saveas(fig6, fullfile(outDIR, plt_nm), 'png');




%% Discard below

%{
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
        

        
        plt_nm = sprintf('compare_io_model_%d.png',idx);
        saveas(fig1, fullfile(saveDIR, plt_nm), 'png');
        
        pause
%}