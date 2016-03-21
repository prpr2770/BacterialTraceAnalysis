function [score_fin, pks, locs] = getChangeScore(sig,leftMargin)
% leftMargin : the index after which you need to measure the MaxValue and
% time-instant of the max-Event. 

alpha = .0;
n = 50; % 2*n+k-2 is the size of the "buffer zone".
k = 10; % window width

try
    score1 = change_detection(sig,n,k,alpha);
    score2 = change_detection(sig(:,end:-1:1),n,k,alpha);
    score2 = score2(end:-1:1); % why do this?
    score_fin = [zeros(1,n-1+k/2), (score1 + score2), zeros(1,n-1+k/2)];
catch
    warning('Error computing change score ')
end


% obtain the maxValue and time instant
% leftMargin = 100; 
T = size(sig,2);
rightMargin = T-leftMargin;
[maxVal, maxIdx] = max(score_fin(1,leftMargin:rightMargin));
%obtain peak values
try
[pks, locs] = findpeaks(score_fin,'MinPeakHeight',0.7*maxVal,'MinPeakProminence',0.1*maxVal);
catch
   warning('Error finding peaks') 
end

subplot(2,1,1)
plot(sig, 'b-', 'linewidth',1);
gridxy(locs,'Color','r','Linestyle',':') ;
subplot(2,1,2)
plot(score_fin, 'r-', 'linewidth',2);
hold on
plot(locs,pks,'*r')
hold off


end

