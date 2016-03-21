function score_fin = getChangeDetectionScore(sig)
% leftMargin : the index after which you need to measure the MaxValue and
% time-instant of the max-Event. 

alpha = .0;
n = 50; % 2*n+k-2 is the size of the "buffer zone".
k = 10; % window width

sprintf('Using change_detection parameters: alpha = %f, n = %d, k = %d ',alpha,n,k);

try
    score1 = change_detection(sig,n,k,alpha);
    score2 = change_detection(sig(:,end:-1:1),n,k,alpha);
    score2 = score2(end:-1:1); % why do this?
    score_fin = [zeros(1,n-1+k/2), (score1 + score2), zeros(1,n-1+k/2)];
catch
    warning('Error computing change score ')
end

end

