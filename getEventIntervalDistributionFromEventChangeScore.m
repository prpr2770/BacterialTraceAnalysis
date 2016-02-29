function allPeakIntervals = getEventIntervalDistributionFromEventChangeScore(changeScores,Fs)

% from collection of ChangeEventScores for each time-trace, extract
% event-interval distribution.
allPeakIntervals = [];
numTraces = size(changeScores,1);

for k = 1:numTraces
data = changeScores(k,:);
[pks, locs] = findpeaks(changeScores,'MinPeakDistance',5);
peakInterval = diff(locs);
allPeakIntervals = [allPeakIntervals (1/Fs*peakInterval)];
end

figure;
histogram(allPeakIntervals);


end