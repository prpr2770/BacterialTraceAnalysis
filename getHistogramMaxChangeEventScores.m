% extract histogram information from changeEventScores

close all;clear all;

try
tracesDirName = 'H:\KraljLab\chgScore\';
tracesFileName = strcat(tracesDirName,'changeScore_Analysis.mat');
load(tracesFileName); 
Fs = 5; % samplingFreq 5Hz
catch
 warning('Error Reading Input Data');
end
%--------------------------------------------------------------------------

N = size(cell_allTraces,2);

EventChangeTimes = [];
maxEventChangeTimes = [];
evtIntervals = [];

data = struct('trace',{},'chgScore',{},'peakLocs',{});

for idx = 1:N
   V = cell_allTraces{idx}; 
   CS = cell_allTraceScores{idx}; 
   T = size(V,2);
   
   % ----------------------------------------------------------------------
   % collate times of peaks
   
   % determine maxVal of changeScore
   leftMargin = 100; rightMargin = T-leftMargin;
   [maxVal, maxInd] = max(CS(1,leftMargin:rightMargin));
   minThreshold = 0.7*maxVal;
   minProminence = 0.1*maxVal;
   
   [pks, locs] = findpeaks(CS,'MinPeakHeight', minThreshold,'MinPeakProminence', minProminence);
   
   % ----------------------------------------------------------------------
   % collate times of peaks

   % extract peak info only in interested region
   indxs = find(locs.*(locs>leftMargin).*(locs<=rightMargin));
   locs = locs(indxs);
   pks = pks(indxs);

   intervals = diff(locs);
   evtIntervals  = [evtIntervals  intervals];
   
   EventChangeTimes = [EventChangeTimes locs]; 
   
   % determine maxVal in the subset region of focus
   [maxVal, maxInd] = max(pks);
   maxEventChangeTimes = [maxEventChangeTimes locs(maxInd)];


   % ARCHIVE VALUES
   % data = struct('trace',{},'chgScore',{},'peakLocs',{});
   data(idx).trace = V;
   data(idx).chgScore = CS;
   data(idx).peakLocs = locs;
end

%archive data-file
save('getHistogramMaxChangeEventScores.mat', 'data')


EventChangeTimes = 1/Fs*EventChangeTimes;
maxEventChangeTimes = 1/Fs*maxEventChangeTimes;

figure;
subplot(1,2,1)
histogram(maxEventChangeTimes, 64)
xlabel('maxEventChange Histogram')
size(maxEventChangeTimes )

subplot(1,2,2)
histogram(EventChangeTimes, 64)
xlabel('EventChange Histogram')

figure;
histogram(evtIntervals,64)