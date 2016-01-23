% =========================================================================
%{
1. Read NormSpectg Batchwise.
2. Compute CodeWords

	Script that does the following:
	+ Randomly select 10 songs files's NORM_SPECTG
	+ Implement k-means, to detect codeWords
	+ Archive CodeWord vectors into a GlobalRepository, for each collection of 10 songs


% numCodeWordsPerBatch = 512;
% songsPerBatch = 10;
% iter = 20;              % kMeans algorithmtracesDirName = 'H:\KraljLab\';


%}
close all;clear all;


% Parameters
numCodeWordsPerBatch = 512; 
songsPerBatch = 200;
iter = 20;                      % for knn-ma-toolbox.

% -------------------------------------------------------------------------
% Folders

tracesDirName = 'H:\KraljLab\';
aggregateDataDirName = fullfile(tracesDirName,'aggSpectG_data\');

% -------------------------------------------------------------------------
% Files

% create file to store Aggregate_Info.
aggData_fileName = strcat(tracesDirName,'aggSPECTG.mat');
aggData_mFile = matfile(aggData_fileName);

% create file to store CodeWords - of BATCHES
cw_batch_fileName = strcat(tracesDirName,'cw_batch_SPECTG.mat');
cw_batch_mFile = matfile(cw_batch_fileName,'Writable',true);

% create file to store CodeWords
cw_fileName = strcat(tracesDirName,'cw_final_SPECTG.mat');
cw_mFile = matfile(cw_fileName,'Writable',true);


% -------------------------------------------------------------------------
% Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
DIMS_SPECTG = aggData_mFile.DIMS_SPECTG ;


% -------------------------------------------------------------------------
% Extract CodeWords




Files = dir(aggregateDataDirName);
rndOrder_songs = randperm(length(Files));

totalBatches = ceil(length(Files)/songsPerBatch);

fileCount = 1;
batchCount = 1;
trcCount = 0;


while (fileCount < length(Files) + 1)
    % ------------------------------------------------------------------------------
    % concatanate all norm_spectg frames for each batch.
    
    BATCH_NORM_SPECTG = [];
    
    % concatanate
    while(mod(fileCount, songsPerBatch) ~=0 && fileCount <= length(Files) )
        warning('Concatanating TRACE data: NORM_SPECTG')
        trcID = rndOrder_songs(fileCount);
        
        if (Files(trcID).bytes > 0 )
            trcCount = trcCount + 1;
            % Read the mat file containing mfcc,dyn_mfcc data
            matfileName=Files(trcID).name;
            spectg_fileName = strcat(aggregateDataDirName,matfileName);
            spectg_mFile = matfile(spectg_fileName);
            
            
            % extract the data needed.
            NORM_SPECTG = spectg_mFile.NORM_SPECTG;
            size(NORM_SPECTG )
            
            % concatanate the frames for the 10 songs
            BATCH_NORM_SPECTG = [BATCH_NORM_SPECTG NORM_SPECTG];
            
        end
        % increment fileCount
        fileCount = fileCount+1
    end
    
    % ------------------------------------------------------------------------------
    % Implement the kmeans algorithm to detect cluster
    size(BATCH_NORM_SPECTG)
    tic
    BATCH_CODEWORDS = ma_kmeans(BATCH_NORM_SPECTG', iter, numCodeWordsPerBatch);
    toc
    BATCH_CODEWORDS = BATCH_CODEWORDS';           % COLVEC - Convert!
    
    % ------------------------------------------------------------------------------
    % store/archive the clusters computed.
    % -------------------------------------------------------
    % archive the norm_dyn_mfcc_ofAllSongs!
    if batchCount > 1
        [nrows ncols] = size(cw_batch_mFile,'CODEWORDS_ALL');
        numFrames = size(BATCH_CODEWORDS,2);
        cw_batch_mFile.CODEWORDS_ALL(:,ncols+1:ncols+numFrames) = BATCH_CODEWORDS;
    else
        cw_batch_mFile.CODEWORDS_ALL = BATCH_CODEWORDS;
    end
    
    
    
    % increment fileCount
    fileCount = fileCount+1
    batchCount = batchCount + 1 
    
    
end

CODEWORDS_ALL = cw_batch_mFile.CODEWORDS_ALL;

FINAL_CODEWORDS = ma_kmeans(CODEWORDS_ALL', iter, numCodeWordsPerBatch);

cw_mFile.CODEWORDS = FINAL_CODEWORDS'; % convert COL-VEC