% Script to do the following:
% 1. Iterate over the Ecoli/Styph/Mix files
% 2. From each, extract windows and the SpectG
% 3. Store Bag of SpectG from each trace. Archive the file.


%%
close all; clear all;

featureType = 'SPECTOGRAM'; % MFCC

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016';
cd(dAnz)

flist = dir('*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz



% Directory to save processed data
saveDir = [dAnz filesep featureType '_' num2str(date)];
if ~exist(saveDir)
    mkdir(saveDir)
end



%% Determine sampleID from fileNames
condList = {'Ecoli';'Styph';'Mix'};
for fID = 1:nfiles
    for j = 1:length(condList)
        if regexp(flist(fID).name,condList{j})
            sampleID(fID) = j;
            break;
        end
    end
end

%% Iterate over each file to do the following:
% 1. Read the denoised signal.
% 2. extract windows and spectG for each window of trace
% 3. Write a new file with the SpectG

% traces_spectG{}; %structure
for fID = 1:nfiles
    fID
    fname = flist(fID).name
    
    % --------------------------------------------------------------------
    % File names to be saved.
    fDIR_windows = [saveDir filesep fname(1:end-4)];
    if ~exist(fDIR_windows)
        mkdir(fDIR_windows)
    end
    
    % create file to store AggregateInfo from all Cells in File.
    aggFeature_fileName = [fDIR_windows filesep 'aggFEATURE.mat'];
    aggFeature_mFile = matfile(aggFeature_fileName,'Writable',true);
    
    
    % --------------------------------------------------------------------
    % Extract data from each file, and process.
    
    fdata = load(fname);
    % extract data
    fields = fieldnames(fdata);
    if numel(fields) == 1
        fields{1};
        data = fdata(1).(fields{1});
    else
        warning('File has more fields than 1.')
    end
    
    % read the denoised signals
    ncells = length(data);
    
    % collecting Aggregate Information
    totalFrames = 0;
    allTraces_MEAN_FEATURE = [];
    
    % iterate over each trace, extract windows.
    for idx = 1:ncells
        
        % obtain signal
        sig = data(idx).sig_dn_wav;
        
        [FEATURE, meanFEATURE] = getFeaturesFromSignal(sig, fs, featureType);
        numFrames = size(FEATURE,2);
        % ------------------------------------------------------------
        % Archive Data for each trace.
        % ------------------------------------------------------------
        
        % create matfile to store Spectrogram Info
        trcName = sprintf('FEATURE_%d.mat',idx);
        feature_fileName = [fDIR_windows filesep trcName];
        feature_mFile = matfile(feature_fileName,'Writable',true);
        
        % save computed values
        feature_mFile.FEATURE = FEATURE;
        feature_mFile.mean_FEATURE = meanFEATURE;
        feature_mFile.numFramesInTrace = numFrames;
        
        % ------------------------------------------------------------
        % Aggregate Information for each File
        % ------------------------------------------------------------
        
        % compute number of frames computed for allCells
        totalFrames = totalFrames + numFrames;
        
        % each trace-mean is stored as a Row-Vector
        if idx ~=1
            allTraces_MEAN_FEATURE = allTraces_MEAN_FEATURE + meanFEATURE;
        else
            allTraces_MEAN_FEATURE = meanFEATURE;
        end
        
        
        
    end
    
    
    % ------------------------------------------------------------
    % Aggregate Information for each File
    % ------------------------------------------------------------
    
    % compute allTraces_MEAN
    MEAN_FEATURE = (1/ncells)* allTraces_MEAN_FEATURE;
    DIMS_FEATURE = length(MEAN_FEATURE);
    
    % save aggregateInfo.
    aggFeature_mFile.TOTAL_FRAMES = totalFrames;
    aggFeature_mFile.MEAN_FEATURE = MEAN_FEATURE;
    aggFeature_mFile.DIMS_FEATURE = DIMS_FEATURE;
    
    
    % -------------------------------------------------------------
    % Compute the Standard-Deviation
    totalFrames = 0;
    varSum = 0;
    for idx = 1:ncells
        
        % create matfile to read Spectrogram Info for each Cell.
        trcName = sprintf('FEATURE_%d.mat',idx);
        feature_fileName = [fDIR_windows filesep trcName];
        feature_mFile = matfile(feature_fileName,'Writable',false);
        
        % extract FEATURE data from the matfile
        FEATURE = feature_mFile.FEATURE;
        numFrames = size(FEATURE,2);
        totalFrames = totalFrames + numFrames;
        
        % compute the variance
        MEAN_FEATURE_mat = repmat(MEAN_FEATURE,1,numFrames);
        varSum = varSum + sum((FEATURE - MEAN_FEATURE_mat).^2,2);  %rowsum
        
        clear MEAN_FEATURE_mat FEATURE numFrames trcName feature_mFile;
    end
    
    
    % Archive the STD_FEATURE information.
    TOTAL_FRAMES = aggFeature_mFile.TOTAL_FRAMES;
    
    if (TOTAL_FRAMES == totalFrames)
        STD_FEATURE = ((1/totalFrames) * varSum).^(0.5);
        aggFeature_mFile.STD_FEATURE = STD_FEATURE;
    else
        error('Error in TOTAL_FRAME count for all songs!')
    end
    
    
end % EndIterateOverFiles

%%
compute_NORM_FEATURES_PROPS