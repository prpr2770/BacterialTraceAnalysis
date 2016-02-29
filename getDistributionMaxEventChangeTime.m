% function [maxEventChangeTime ,allPeakIntervals] = getDistributionMaxEventChangeTime(traces,Fs, k)
% script that determines histogram of the time-instance of maxChangeEvent
% detection algo, for all traces.


alpha = .0;
n = 50; % 2*n+k-2 is the size of the "buffer zone".
% k = 10; % window width


numTraces = size(traces,1);
maxEventChangeTime = zeros(1,numTraces);
allPeakIntervals = [];

%create structure to store changeEventScore for each trace.

for idx = 1:numTraces
    idx
    sig = traces(idx,:);
    try
        score1 = change_detection(sig,n,k,alpha);
        score2 = change_detection(sig(:,end:-1:1),n,k,alpha);
        score2 = score2(end:-1:1); % why do this?
        score_fin = [zeros(1,n-1+k/2), (score1 + score2), zeros(1,n-1+k/2)];
        size(score_fin)
        
        % obtain the maxValue and time instant
        [maxVal, maxIdx] = max(score_fin);
        maxEventChangeTime(idx) = maxIdx;
        
        %obtain peak values
        [pks, locs] = findpeaks(score_fin,'MinPeakHeight',maxVal/20,'MinPeakProminence',maxVal/10);
        peakInterval = diff(locs);
        
        allPeakIntervals = [allPeakIntervals (1/Fs*peakInterval)];
        
        % archive intermediate results
        cell_allTraces{idx} = sig;
        cell_allTraceScores{idx} = score_fin;
        cell_allPeakIntervals{idx} = (1/Fs*peakInterval);
        cell_allPeakLocs{idx} = locs;
        
        % plots
        fig1 = figure(1)
        subplot(2,1,1)
        plot(sig, 'b-', 'linewidth',1);
        gridxy(locs,'Color','r','Linestyle',':') ;
        subplot(2,1,2)
        plot(score_fin, 'r-', 'linewidth',2);
        hold on
        plot(locs,pks,'*r')
        hold off
        
        tracesDirName = 'H:\KraljLab\chgScore\';
        fn = sprintf('trackID_%d.png',idx);
        plotFileName = strcat(tracesDirName,fn);
        saveas(fig1, fullfile(tracesDirName, fn), 'png');
        close all
    catch
        strError = sprintf('Error Processing Trace: %d',idx);
        warning(strError);
    end
end

maxEventChangeTime = 1/Fs * maxEventChangeTime;

figure;
subplot(1,2,1)
histogram(maxEventChangeTime,32);
xlabel('maxChangeEventTime(s)')

subplot(1,2,2)
histogram(allPeakIntervals,32);
xlabel('event-interval(s)')

set(gcf,'NextPlot','add');
axes;
strTitle = sprintf('Window size: %d',k);
h = title(strTitle);
set(gca,'Visible','off');
set(h,'Visible','on');
savefig('HistogramEventDurations.fig')

save('changeScore_Analysis.mat','cell_allTraces','cell_allTraceScores','cell_allPeakIntervals','cell_allPeakLocs')










% end
