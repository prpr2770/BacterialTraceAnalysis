% To be called after the extraction of timeTraces from Video.

% Compare the time-domain signals for each folder. Also, plot the
% aggregate feature vectors.

% Script to do the following:
% 1. Iterate over the Ecoli/Styph/Mix files
% 2. From each, extract windows and the SpectG
% 3. Store Bag of SpectG from each trace. Archive the file.


%%
close all; clear all;

featureType = 'VISUALIZE';%'SPECTOGRAM'; % MFCC
signalType = 'VOLTAGE';

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016';
cd(dAnz)

flist = dir('*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz



% Directory to save processed data
saveDir = [dAnz filesep featureType '_' num2str(date)];
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
    
    fParams = [signalType '_' condList{sampleID(fID)} '_'];
    
    % --------------------------------------------------------------------
    % File names to be saved.
    fDIR = [saveDir filesep fname(1:end-4)];
    if ~exist(fDIR)
        mkdir(fDIR)
    end
    
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
    
    
    % iterate over each trace, extract windows.
    all_sig = [];
    all_sig_dn = [];
    
    all_sig_chg = [];
    all_sig_dn_chg = [];
    
    for idx = 1:ncells
        
        % obtain signal
        sig_dn = data(idx).sig_dn_wav;
        sig = data(idx).sig;
        
        sig_chg = data(idx).sig_chgScr;
        sig_dn_chg = data(idx).sig_dn_wav_chgScr;
        
        % compile all signals
        all_sig = [all_sig; sig];
        all_sig_dn = [all_sig_dn; sig_dn];
        
        all_sig_chg = [all_sig_chg; sig_chg/max(sig_chg)];
        all_sig_dn_chg = [all_sig_dn_chg; sig_dn_chg/max(sig_dn_chg)];
% % %         %%
% % %         
% % %         % -------------------------------------------------------------------------
% % %         % Fig1 = Subplot(2,1)
% % %         % 1. Volt : noisy -vs- denoisy
% % %         % 2. Volt events  : noisy-vs-dn
% % %         
% % %         fig1 = figure(1)
% % %         subplot(2,1,1)
% % %         p1 = plot(sig,'g'); hold on;
% % %         p2 = plot(sig_dn,'r'); hold off;
% % %         legend([p1,p2],'nsy','de-nsy')
% % %         ylabel('voltage')
% % %         subplot(2,1,2)
% % %         plot(sig_chg/max(sig_chg),'g'); hold on;
% % %         plot(sig_dn_chg/max(sig_dn_chg),'r'); hold off;
% % %         ylabel('chg-score')
% % %         
% % %         set(gcf,'NextPlot','add');
% % %         axes;
% % %         heading = [signalType sprintf('signals for trackID:%d',idx)];
% % %         h = title(heading);
% % %         set(gca,'Visible','off');
% % %         set(h,'Visible','on');
% % %         
% % %         plt_nm = [signalType sprintf('_trackID_%d.png',idx)];
% % %         saveas(fig1, fullfile(fDIR, plt_nm), 'png');
% % %         
% % %         % =========================================================================
% % %         % 3. Iterate over each trace and plot the volt-ca-signals
% % %         
% % %         close all;
        
        
        
    end
    
    mean_sig = mean(all_sig,1);
    mean_sig_dn = mean(all_sig_dn,1);
    
    
    fig2 = figure(2)
    subplot(3,1,1)
    p1 = plot(all_sig','Color',[0.8 0.8 0.8]); hold on;
    p2 = plot(mean_sig,'Color',[0.8 0.4 0.8]); hold off;
    ylabel(signalType)
    subplot(3,1,2)
    p1 = plot(all_sig_dn','Color',[0.8 0.8 0.8]); hold on;
    p2 = plot(mean_sig_dn,'Color',[0.4 0.8 0.8]); hold off;
    ylabel([signalType '-denoised'])
    subplot(3,1,3)
    p1 = plot(mean_sig,'Color',[0.8 0.4 0.8]); hold on;
    p2 = plot(mean_sig_dn,'Color',[0.4 0.8 0.8]); hold off;
    ylabel([signalType ' mean :(de)noisy'])
    
    set(gcf,'NextPlot','add');
    axes;
    heading = [signalType ': aggregate and mean'];
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    plt_nm = [fParams '_aggregate'];
    saveas(fig2, fullfile(fDIR, plt_nm), 'png');
    
    %%
    fig3 = figure(3)
    colormap jet
    subplot(1,2,1)
    imagesc(all_sig)
    subplot(1,2,2)
    imagesc(all_sig_dn)
    colorbar
    
    set(gcf,'NextPlot','add');
    axes;
    heading = [signalType ': aggregate matrix'];
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    
    plt_nm = [fParams '_aggMatrix'];
    saveas(fig3, fullfile(fDIR, plt_nm), 'png');
    
    %%
    fig4 = figure(4)
    colormap jet
    subplot(1,2,1)
    imagesc(all_sig_chg)
    colorbar
    subplot(1,2,2)
    imagesc(all_sig_dn_chg)
    
    set(gcf,'NextPlot','add');
    axes;
    heading = [signalType ': ChgScr aggMatrix'];
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
  
    
      plt_nm = [fParams '_chgScr_aggMatrix'];
    saveas(fig4, fullfile(fDIR, plt_nm), 'png');
    
%%

fig5 = figure(5)
mesh(all_sig_dn_chg)
title('agg-chgScr-denoised')
plt_nm = [fParams '_chgScr_aggMatrix'];
savefig(fig5, fullfile(fDIR, plt_nm))

fig6 = figure(6)
mesh(all_sig_dn)
title('agg-sig-denoised')
plt_nm = [fParams '_sig_aggMatrix'];
savefig(fig6, fullfile(fDIR, plt_nm))

%%
close all;
end % EndIterateOverFiles

%%
