% Script to do the following:
% 1. Iterate over the Ecoli/Styph/Mix files
% 2. From each, extract windows and the SpectG
% 3. Store Bag of SpectG from each trace. Archive the file.


%%
close all; clear all;

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016';
cd(dAnz)

flist = dir('*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz

% Directory to save processed data
saveDir = [dAnz filesep 'spectG_' num2str(date)];
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
    aggSpectg_fileName = [fDIR_windows filesep 'aggSPECTG.mat'];
    aggSpectg_mFile = matfile(aggSpectg_fileName,'Writable',true);
    
    
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
    allTraces_MEAN_SPECTG = [];
    
    % iterate over each trace, extract windows.
    for idx = 1:ncells
        
        % obtain signal
        sig = data(idx).sig_dn_wav;
        
        % window details
        Nx = length(sig);
        nsc = floor(Nx/4.5);
        novlap = floor(3*nsc/4);
        nfft = max(256,2^nextpow2(nsc));
        wind = hamming(nsc);
        
        % ------------------------------------------------------------
        % Compute Spectorgram
        % ------------------------------------------------------------
        
        % compute Spectrogram
        [SPECTG, F, T] = spectrogram(sig,wind,novlap,nfft,fs); %spectrogram(sig,wind,novlap,nfft,fs);
        SPECTG = log10(abs(SPECTG)); % take the log10 of the magnitude of Spectrogram.
        
        % compute the MeanSpec of song.
        numFrames = size(SPECTG,2);
        meanSPECTG = (1/numFrames)*sum(SPECTG,2);
        
        % ------------------------------------------------------------
        % Archive Data for each trace.
        % ------------------------------------------------------------
        
        % create matfile to store Spectrogram Info
        trcName = sprintf('SPECTG_%d.mat',idx);
        spectg_fileName = [fDIR_windows filesep trcName];
        spectg_mFile = matfile(spectg_fileName,'Writable',true);
        
        % save computed values
        spectg_mFile.SPECTG = SPECTG;
        spectg_mFile.F = F;
        spectg_mFile.T = T;
        spectg_mFile.mean_SPECTG = meanSPECTG;
        spectg_mFile.numFramesInTrace = numFrames;
        
        % ------------------------------------------------------------
        % Aggregate Information for each File
        % ------------------------------------------------------------
        
        % compute number of frames computed for allCells
        totalFrames = totalFrames + numFrames;
        
        % each trace-mean is stored as a Row-Vector
        if idx ~=1
            allTraces_MEAN_SPECTG = allTraces_MEAN_SPECTG + meanSPECTG;
        else
            allTraces_MEAN_SPECTG = meanSPECTG;
        end
        
        
        
    end
    
    
    % ------------------------------------------------------------
    % Aggregate Information for each File
    % ------------------------------------------------------------
    
    % compute allTraces_MEAN
    MEAN_SPECTG = (1/ncells)* allTraces_MEAN_SPECTG;
    DIMS_SPECTG = length(MEAN_SPECTG);
    
    % save aggregateInfo.
    aggSpectg_mFile.TOTAL_FRAMES = totalFrames;
    aggSpectg_mFile.MEAN_SPECTG = MEAN_SPECTG;
    aggSpectg_mFile.DIMS_SPECTG = DIMS_SPECTG;
    
    
    % -------------------------------------------------------------
    % Compute the Standard-Deviation
    totalFrames = 0;
    varSum = 0;
    for idx = 1:ncells
        
        % create matfile to read Spectrogram Info for each Cell.
        trcName = sprintf('SPECTG_%d.mat',idx);
        spectg_fileName = [fDIR_windows filesep trcName];
        spectg_mFile = matfile(spectg_fileName,'Writable',false);
        
        % extract SPECTG data from the matfile
        SPECTG = spectg_mFile.SPECTG;
        numFrames = size(SPECTG,2);
        totalFrames = totalFrames + numFrames;
        
        % compute the variance
        MEAN_SPECTG_mat = repmat(MEAN_SPECTG,1,numFrames);
        varSum = varSum + sum((SPECTG - MEAN_SPECTG_mat).^2,2);  %rowsum
        
        clear MEAN_SPECTG_mat SPECTG numFrames trcName spectg_mFile;
    end
    
    
    % Archive the STD_SPECTG information.
    TOTAL_FRAMES = aggSpectg_mFile.TOTAL_FRAMES;
    
    if (TOTAL_FRAMES == totalFrames)
        STD_SPECTG = ((1/totalFrames) * varSum).^(0.5);
        aggSpectg_mFile.STD_SPECTG = STD_SPECTG;
    else
        error('Error in TOTAL_FRAME count for all songs!')
    end
    
    
end % EndIterateOverFiles

%%
compute_Norm_SpectG_PROPS_mixture