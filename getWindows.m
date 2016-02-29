function [win_sigs, winCount] = getWindows(sigs, W, delay)
% function to extract windows from time-series.
% input: sigs - row-vector time-series.
% W: window length
% delay : delay between moving windows. 
% 
% Output:
% win_sigs : col-vectors of windows.
% winCounts : totalNumber of windows extracted. 
% -------------------------------------------------------------------

numSamples = size(sigs,2);
numTraces = size(sigs,1);

% sample traces as col-vectors
trace = sigs';

lastWindowStart = numSamples - W + 1;

idx = 1;
winCount = 0;

win_trace = [];
% win_fft_trace = [];


while idx <=lastWindowStart
    winCount = winCount + 1;
    window = trace(idx:idx+W-1,:);
    
    % extract windowed info.
    win_trace(:,:,winCount) = window;
%     win_fft_trace(:,:,winCount) = fft(window,W);
    
    % increment index
    idx = idx + delay;
end

win_sigs = win_trace;
end


%  Unit-test
% % sigs = rand(2,10); W = 5 ; delay = 2;
% % wind_sigs = getWindows(sigs, W, delay)