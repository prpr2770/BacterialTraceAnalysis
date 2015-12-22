% =========================================================================
%{
1. Create SPECTG_trcID.mat and agg_SPECTG.mat
2. In allSPECTG.mat store:
    a) SPECTG       : for each song. 
    b) mean_SPECTG  : for each song. 
    a) NORM_SPECTG  : for each song
3. agg_SPECTG.mat
    a) MEAN_SPECTG  : for all frames!
    b) STD_SPECTG   : for all frames!
    c) TOTAL_FRAMES :

%}
close all;clear all;


tracesDirName = 'H:\KraljLab\';
tracesFileName = strcat(tracesDirName,'PROPS_data.mat');
load(tracesFileName);
fs = 5; % samplingFreq 5Hz

% ----------------------------------------------------------------------------------------------

% data directory to store spectrogram info. 
aggregateDataDirName = fullfile(tracesDirName,'aggSpectG_data\');
status = mkdir(aggregateDataDirName)

% create file to store Norm_Spectrogram
aggSpectg_fileName = strcat(tracesDirName,'aggSPECTG.mat');
aggSpectg_mFile = matfile(aggSpectg_fileName,'Writable',true);

% ----------------------------------------------------------------------------------------------

[numTraces, numSamples] = size(intens);

totalFrames = 0;
allTraces_MEAN_SPECTG = [];
for trcID = 1:numTraces
    sig = intens(trcID,:)';

    Nx = length(sig);
    nsc = floor(Nx/4.5);
    novlap = floor(3*nsc/4);
    nfft = max(256,2^nextpow2(nsc));
    wind = hamming(nsc);
    
    % compute Spectrogram
    [SPECTG, F, T] = spectrogram(sig,wind,novlap,nfft,fs);
    SPECTG = log10(abs(SPECTG)); % take the log10 of the magnitude of Spectrogram.
    
    % compute the MeanSpec of song.
    numFrames = size(SPECTG,2);
    meanSPECTG = (1/numFrames)*sum(SPECTG,2);
    
    % write all the data
    totalFrames = totalFrames + numFrames;

    % create matfile to store Spectrogram Info
    trcName = sprintf('SPECTG_%d.mat',trcID);
    spectg_fileName = strcat(aggregateDataDirName,trcName);
    spectg_mFile = matfile(spectg_fileName,'Writable',true);
    
    % save computed values
    spectg_mFile.SPECTG = SPECTG;
    spectg_mFile.F = F;
    spectg_mFile.T = T;
    spectg_mFile.mean_SPECTG = meanSPECTG;
    spectg_mFile.numFramesInTrace = numFrames;
    
    % each trace-mean is stored as a Row-Vector
    if trcID ~=1
        allTraces_MEAN_SPECTG = allTraces_MEAN_SPECTG + meanSPECTG;
    else
        allTraces_MEAN_SPECTG = meanSPECTG; 
    end
end

% compute allTraces_MEAN
MEAN_SPECTG = (1/numTraces)* allTraces_MEAN_SPECTG;
DIMS_SPECTG = length(MEAN_SPECTG);

% save aggregateInfo.
aggSpectg_mFile.TOTAL_FRAMES = totalFrames;
aggSpectg_mFile.MEAN_SPECTG = MEAN_SPECTG;
aggSpectg_mFile.DIMS_SPECTG = DIMS_SPECTG;
% ----------------------------------------------------------------------------------------------
% Compute STD_SPECTG:

varSum = zeros(DIMS_SPECTG,1);
totalFrames = 0;

trcID = 0;
Files = dir(aggregateDataDirName);


for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
    if (Files(k).bytes >0 )
        trcID = trcID + 1;
        
        % Read the mat file containing SPECTG for the TRACE
        trace_fileName = Files(k).name;
        trace_fileName = strcat(aggregateDataDirName,trace_fileName);
        trace_mFile = matfile(trace_fileName);
        

        % extract the data needed.
        SPECTG = trace_mFile.SPECTG;
        [coeffs numFrames] = size(SPECTG);
        totalFrames = totalFrames + numFrames;

        % compute the variance
        MEAN_SPECTG_mat = repmat(MEAN_SPECTG,1,numFrames);
        varSum = varSum + sum((SPECTG - MEAN_SPECTG_mat).^2,2);  %rowsum

        clear MEAN_SPECTG_mat SPECTG;
        
    end
end
TOTAL_FRAMES = aggSpectg_mFile.TOTAL_FRAMES;

if (TOTAL_FRAMES == totalFrames)
    STD_SPECTG = ((1/totalFrames) * varSum).^(0.5);
    aggSpectg_mFile.STD_SPECTG = STD_SPECTG;
else
    error('Error in TOTAL_FRAME count for all songs!')
end

% ----------------------------------------------------------------------------------------------
% Compute norm_SPECTG:

trcID = 0;
Files = dir(aggregateDataDirName);

for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
    if (Files(k).bytes >0 )
        trcID = trcID + 1;
        
        % Read the mat file containing SPECTG for the TRACE
        trace_fileName = Files(k).name;
        trace_fileName = strcat(aggregateDataDirName,trace_fileName);
        trace_mFile = matfile(trace_fileName,'Writable',true);
        

        % extract the data needed.
        SPECTG = trace_mFile.SPECTG;
        [coeffs numFrames] = size(SPECTG);
        totalFrames = totalFrames + numFrames;

        % compute the variance
        MEAN_SPECTG_mat = repmat(MEAN_SPECTG,1,numFrames);
        STD_SPECTG_mat = repmat(STD_SPECTG,1,numFrames);

        % compute Normalized SPECTG
        NORM_SPECTG = (SPECTG - MEAN_SPECTG_mat)./STD_SPECTG_mat;
        
        % archive computed value.
        trace_mFile.NORM_SPECTG = NORM_SPECTG;
        
        % remove 
        clear NORM_SPECTG STD_SPECTG_mat MEAN_SPECTG_mat SPECTG;
        
    end
end

