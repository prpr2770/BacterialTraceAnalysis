function [subseq, subseq_tseries, fh] = getSubsequence(seq,tIndx,leftBorder,rightBorder)
% function to extract a subsequence centered at a location given in tIndx,
% and ranging from (indx-leftMargin, indx+rightMargin) of the seq.

seqT = length(seq);
subseq = {};
interval_times = [];
for k = 1:length(tIndx)
    t = tIndx(k);
    ti = t-leftBorder;
    tf = t + rightBorder;
    
    if ti < 1 %ti = ti*(ti>=1) + 1*(ti<1);
        ti = 1;
    elseif tf > seqT %ti = tf*(tf<=seqT) + seqT*(tf>seqT);
        tf = seqT;
    end
    
    subseq{k}.sig = seq(ti:tf);
    subseq{k}.time = ti:tf;
    
    interval_times = [interval_times ti:tf];

end

% plotting chosen intervals of the sequence
interval_times = unique(sort(interval_times));
Indicator = zeros(size(seq));
Indicator(interval_times) = ones(size(interval_times));
subseq_tseries = seq.*Indicator;

% % % ----------------------------------
% % % plot 
% % % ----------------------------------
% % fh = figure;
% % plot(seq,'g')
% % hold on
% % plot(subseq_tseries,'r')
% % hold off
% % xlabel('time')
% % % ----------------------------------

end

