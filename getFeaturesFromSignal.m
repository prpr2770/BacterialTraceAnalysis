function [FEATURE, meanFEATURE] = getFeaturesFromSignal(sig, fs, featureType)
% Function that extracts different features from a signal, based on
% featureType


switch featureType
    case 'SPECTROGRAM'
    FEATURE= spectrogram(sig,fs); %spectrogram(sig,wind,novlap,nfft,fs);
    FEATURE = log10(abs(FEATURE)); % take the log10 of the magnitude of Spectrogram.
    % compute the MeanSpec of song.
        numFrames = size(FEATURE,2);
        meanFEATURE = (1/numFrames)*sum(FEATURE,2);
        
	case 'FFT'
	FEATURE = abs(fft(sig))';
	meanFEATURE = FEATURE;
	
		
    case 'WIN_SPECTROGRAM'
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
      
        p = struct('fs',fs,'visu',0,'fft_size',128,'hopsize',64,'num_ceps_coeffs',14,'mel_filt_bank',[0.1 2.5 10],'dB_max',24);

        % generate mfcc and dct
        [FEATURE, DCT] = ma_mfcc_bio(sig',p);
        numFrames = size(FEATURE,2);        %each frame is a column-vector
        meanFEATURE = (1/numFrames)*sum(FEATURE,2);

      
    
        
    otherwise
        error('NO VALID FEATURE TYPE SELECTED.')
end

end