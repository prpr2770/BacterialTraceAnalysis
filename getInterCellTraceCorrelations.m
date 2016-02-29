% extract information regarding the correlation between the volt-volt and
% ca-ca signals for all cells.

clear all;
close all;
% clear all;

loadData = 1;
if loadData
    try
        inDirName = 'H:\KraljLab\voltCal_results\';
        fName = strcat(inDirName ,'changeScore_Analysis.mat');
        load(fName);
        totalsigs = size(data,2);
        loadData = 0;
    catch
        error('Data loading error.')
    end
end

outDirName = 'H:\KraljLab\voltCal_results\intercellXcorr\';

count = 0;

all_V = [];
all_C = [];
all_chgV = [];
all_chgC = [];


for idx = 1:totalsigs
    idx
    chgV = data(idx).voltChangeScore;
    chgC = data(idx).calcChangeScore;
    V = data(idx).TimeTraceVolt;
    C = data(idx).TimeTraceCa;
    
    % % % % +==============================================================
    % extract MEAN and then Normalize.
    V = getNormalized(V);
    C = getNormalized(C);
    
    chgV = getNormalized(chgV);
    chgC = getNormalized(chgC);
    
    
    if(idx>1 && (length(all_V) == length(V)) &&  (length(all_C) == length(C)))
        all_V = [all_V; V];
        all_C = [all_C; C];
        
        all_chgV = [all_chgV; chgV];
        all_chgC = [all_chgC; chgC];
        
        count = count + 1;
        
    elseif idx ==1
        all_V = [all_V; V];
        all_C = [all_C; C];
        
        
        all_chgV = [all_chgV; chgV];
        all_chgC = [all_chgC; chgC];
        
        count = count + 1;
    end
    
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Compute the cross-correlations.

N = count;

for i = 1:N
    v1 = all_V(i,:);
    c1 = all_C(i,:);
    
    xcorr_v12_sum = [];
    xcorr_c12_sum = [];
    
    tau_volt = zeros(1,count);
    tau_ca = zeros(1,count);
    
    for j=1:N
        v2 = all_V(j,:);
        c2 = all_C(j,:);
        
        [v12_corr, v12_lag] = xcorr(v1,v2);
        [~, tID] = max(v12_corr);
        tau_volt(j) = v12_lag(tID);
        
        [c12_corr, c12_lag] = xcorr(c1,c2);
        [~, tID] = max(c12_corr);
        tau_ca(j) = c12_lag(tID);
        
        
        xcorr_v12_sum = [xcorr_v12_sum; v12_corr];
        xcorr_c12_sum = [xcorr_c12_sum; c12_corr];
        
        
%         ----------------------------
%         PLOTS
%         ----------------------------
        fig1 = figure(1)
        subplot(2,2,1)
        plot(v1,'r')
        hold on
        plot(v2,'g')
        hold off
        ylabel('volt')
        
        subplot(2,2,3)
        plot(c1,'r')
        hold on
        plot(c2,'g')
        hold off
        ylabel('calcium')
        
        
        subplot(2,2,2)
        plot(v12_lag, v12_corr )
        title('xcorr V1-v2')
        
        subplot(2,2,4)
        plot(c12_lag, c12_corr )
        title('xcorr C1-C2')
%         ------------------------------
%         save plot figure
       
        fn = sprintf('trackID_%d_%d.png',i,j);
        saveas(fig1, fullfile(outDirName, fn), 'png');
        close all
        
        
    end
    
    xcorr_v12_mean = sum(xcorr_v12_sum)/count;
    xcorr_c12_mean = sum(xcorr_c12_sum)/count;
    
    % Plot figures for given v1-c1 combination
    fig1 = figure(1)
    subplot(2,1,1)
    plot(v12_lag, xcorr_v12_mean)
    xlabel('avg. Volt xcorr')
    subplot(2,1,2)
    plot(c12_lag, xcorr_c12_mean)
    xlabel('avg Ca xcorr')
    fn = sprintf('avg_xcorr_trackID_%d.png',i);
    saveas(fig1, fullfile(outDirName, fn), 'png');
    
    
    fig2 = figure(2);
    subplot(2,1,1)
    plot(tau_volt)
    ylabel('Volt(\tau)')
    subplot(2,1,2)
    plot(tau_ca)
    ylabel('Ca(\tau)')
    title('\tau-delay from cross-correlation')
    fn = sprintf('tauDelay_trackID_%d.png',i);
    saveas(fig2, fullfile(outDirName, fn), 'png');
    
    % plot the correlations for all signals as an image
    fig3 = figure(3);
    subplot(2,1,1)
    imagesc(xcorr_v12_sum)
    colorbar
    title('xcorr Volt')
    subplot(2,1,2)
    imagesc(xcorr_c12_sum)
    colorbar
    title('xcorr Ca')
    fn = sprintf('all_xcorr_trackID_%d.png',i);
    saveas(fig3, fullfile(outDirName, fn), 'png');
    
    close all;
    
    % Archive values
    xcorrData(i).voltTrace = v1;
    xcorrData(i).caTrace = c1;
    xcorrData(i).xcorrAllVolt = xcorr_v12_sum;
    xcorrData(i).xcorrAllCa = xcorr_c12_sum;
    xcorrData(i).voltLag = v12_lag;
    xcorrData(i).caLag = c12_lag;
end

% save the xcorrData file!
fname = 'xcorrData.mat';
save(fullfile(outDirName, fname),'xcorrData');



