% Script to do the following:
% 1. Iterate over the Ecoli/Styph/Mix files
% 2. from each file, read the dn_sig and dn_sig_chgScr. 
% 3. detect location of max chgScr, and extract motif. 


%%
close all; clear all;

scriptTask = 'MOTIFS';%'SPECTOGRAM'; % MFCC
signalType = 'CALCIUM';

dAnz = 'H:\KraljLab\2016-02-18-PROPS_CALC_Ecoli_vs_Salmonella\ChangeScores';
cd(dAnz)

flist = dir('*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz



% Directory to save processed data
saveDir = [dAnz filesep scriptTask '_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end



%% Determine sampleID from fileNames
condList = {'Ecoli';'Styph';'Mix'};
for fID = 1:nfiles
    for j = 1:length(condList)
        if regexp(flist(fID).name,condList{j})
            sampleID(fID) = j;
            break;
        end
    end
end

%% Iterate over each file to do the following:
% 1. Read the denoised signal.
% 2. extract windows and spectG for each window of trace
% 3. Write a new file with the SpectG

% traces_spectG{}; %structure
for fID = 1:nfiles
    fID
    fname = flist(fID).name
    
    % --------------------------------------------------------------------
    % File names to be saved.
    fDIR = [saveDir filesep fname(1:end-4)];
    if ~exist(fDIR)
        mkdir(fDIR)
    end
    
    % create file to store AggregateInfo from all Cells in File.
    motifs_fileName = [fDIR filesep 'motifs.mat'];
    motifs_mFile = matfile(motifs_fileName,'Writable',true);
    
    
    % --------------------------------------------------------------------
    % Extract data from each file, and process.
    
    fdata = load(fname);
    % extract data
    fields = fieldnames(fdata);
    if numel(fields) == 1
        fields{1};
        data = fdata(1).(fields{1});
    else
        warning('File has more fields than 1.')
    end
    
    % read the denoised signals
    ncells = length(data);
    
    all_motifs = [];
    % iterate over each trace, extract windows.
    for idx = 1:ncells
        
        % obtain signal
        sig = data(idx).sig_dn_wav;
        sig_chg = data(idx).sig_dn_wav_chgScr;
        
        % determine max sig_chg location
        [maxVal, maxIdx] = max(sig_chg);
        
        % extract motif
        sig_motif = sig(maxIdx-20:maxIdx+40);
        
        znrm_motif = (sig_motif - mean(sig_motif))/var(sig_motif);
        
        % compile list of motifs
        all_motifs = [all_motifs; znrm_motif];
    
    end
    
    % extract eigen-motifs
    source = all_motifs;
    target = source;
    k = 32;
    [S, W]= GeneralizedGreedySelection(target, source, k)
    eig_motifs = source(S,:);
    
    % archive the motifs and eigMotifs
        motifs_mFile.ALL_MOTIFS = all_motifs;
        motifs_mFile.EIG_MOTIFS = eig_motifs;
        
        %%
        % visualize motifs
        
        
    fig1 = figure(1)
    numPlots = length(S);
    numRows = min(4, length(S));
    numCols = ceil(numPlots/numRows);
    
    for idy =1:length(S)
        subplot(numRows,numCols,idy)
        p1 = plot(eig_motifs(idy,:),'b'); hold on;
        
    end
    
    
    set(gcf,'NextPlot','add');
    axes;
    heading = [signalType 'Motifs: from ' condList{sampleID(fID)}];
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    plt_nm = ['eigMotifs_' condList{sampleID(fID)}];
    saveas(fig1, fullfile(fDIR, plt_nm), 'png');
    
    %%
    fig2 = figure(2)
    plot(all_motifs','Color',[0.8 0.8 0.8]);
    hold on
    plot(mean(all_motifs,1),'Color',[0.8 0.4 0.8]);
    hold off
    plt_nm = ['allMotifs_' condList{sampleID(fID)}];
    saveas(fig2, fullfile(fDIR, plt_nm), 'png');
    
    
    %%
        close all;
        
end % EndIterateOverFiles

