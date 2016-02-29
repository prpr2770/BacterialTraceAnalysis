function features_win_sigs = getFeatureRepresentation(win_sigs,winLen)
% script: extract feature-vector representation of windows from
% time-traces. 

numWindows = size(win_sigs,3);
features_win_sigs = zeros(size(win_sigs));

for idx = 1: numWindows
    % temp: fft of windowed signal.
    features_win_sigs(:,:,idx) = abs(fft(win_sigs(:,:,idx),winLen));
    
    % alternative: MFCC of signal!
end


end