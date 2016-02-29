function [Traces, delays, x1_delayed, x2_delayed, corr_x1_dx1, corr_x1_dx2, corr_x2_dx1, corr_x2_dx2] = getCorr_OfTracePair(traces,Fs)
% normalize the traces, compute the correlation and auto-correlation.
% plot the same

traces;
Ts = 1/Fs;

% at zero-delay
normTraces = getNormalizedData_ColVecs(traces);
autoCorr = sum(normTraces.*normTraces);
xCorr = sum(normTraces(:,1).*normTraces(:,2));
Corr = [autoCorr xCorr];

% ------------------------------------------------------------------------
% incorporating Delays

totalSamples = length(traces);
totalTraceLength = Ts*totalSamples;

% Assume 1/20th of totalTrace at the START and END can be overlooked.
% We shall analyze correlation only in-between.
marginLength = floor(totalSamples/20)
dataLength = totalTraceLength - 2*marginLength;

% divide margin into atleast 5 units
if marginLength > 10
    shiftLength = floor(marginLength/10);
else
    shiftLength = 1;
    error('Error in Margin and interval/delay creation.')
end


x1_delayed = [];
x2_delayed = [];

delays = (-8:1:8);
numDelayedVecs = length(delays);
for lag = 1:length(delays)
    k = delays(lag);
    indx = marginLength + k*shiftLength;
    x1_delayed = [x1_delayed traces((indx+1):(indx+dataLength),1)];
    x2_delayed = [x2_delayed traces((indx+1):(indx+dataLength),2)];
    
    % identify 0-shifted signal
    if k==0
        x1 = traces((indx+1):(indx+dataLength),1);
        x2 = traces((indx+1):(indx+dataLength),2);
    end
end

% identify normalized vectors whose correlations are to be computed
Traces = [x1 x2];
normTraces = getNormalizedData_ColVecs(Traces);
norm_x1_delayed = getNormalizedData_ColVecs(x1_delayed);
norm_x2_delayed = getNormalizedData_ColVecs(x2_delayed);

% computing correlations
norm_x1 = normTraces(:,1);
corr_x1_dx1 = sum(repmat(norm_x1,1,numDelayedVecs).*norm_x1_delayed);
corr_x1_dx2 = sum(repmat(norm_x1,1,numDelayedVecs).*norm_x2_delayed);

norm_x2 = normTraces(:,2);
corr_x2_dx1 = sum(repmat(norm_x2,1,numDelayedVecs).*norm_x1_delayed);
corr_x2_dx2 = sum(repmat(norm_x2,1,numDelayedVecs).*norm_x2_delayed);

% plot the correlation functions
figure;
subplot(2,2,1)
plot(delays,corr_x1_dx1)
xlabel('lag interval')
title('autocorr: Trace 1')
subplot(2,2,2)
plot(delays,corr_x2_dx2)
title('autocorr: Trace 2')
xlabel('lag interval')

% cross-correlations
subplot(2,2,3)
plot(delays,corr_x1_dx2)
title('xcorr: Trc1 + delayed-Trc2')
xlabel('lag interval')

subplot(2,2,4)
plot(delays,corr_x2_dx1)
title('xcorr: Trc2 + delayed-Trc1')
xlabel('lag interval')

end