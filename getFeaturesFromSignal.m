function [FEATURE, meanFEATURE] = getFeaturesFromSignal(sig, fs, featureType)
% Function that extracts different features from a signal, based on
% featureType


switch featureType
    case 'SPECTOGRAM'
    %% window details
        Nx = length(sig);
        nsc = floor(Nx/4.5);
        novlap = floor(3*nsc/4);
        nfft = max(256,2^nextpow2(nsc));
        wind = hamming(nsc);
        
        % compute Spectrogram
        [FEATURE, F, T] = spectrogram(sig,wind,novlap,nfft,fs); %spectrogram(sig,wind,novlap,nfft,fs);
        FEATURE = log10(abs(FEATURE)); % take the log10 of the magnitude of Spectrogram.
        
        % compute the MeanSpec of song.
        numFrames = size(FEATURE,2);
        meanFEATURE = (1/numFrames)*sum(FEATURE,2);
      
    case 'MFCC'
      %%
      
      
        
    otherwise
        error('NO VALID FEATURE TYPE SELECTED.')
end

end