% =========================================================================
%{
1. Compute CodeWords_Histogram for each trace.
2. Generate CodeWord_Histogram File -> containing all Songs.

Script that does the following:
+ Reads all the CODEWORDS
+ Reads every individual song, and obtain the cw-historgram representation of each song.
+ Save the code-words of each song in a different-file.


%}
close all;clear all;


% Parameters
iter = 20;                      % for knn-ma-toolbox.
tau = 2;
numCodeWordsPerBatch = 512; 

% -------------------------------------------------------------------------
% Folders

tracesDirName = 'H:\KraljLab\';
aggregateDataDirName = fullfile(tracesDirName,'aggSpectG_data\');

% -------------------------------------------------------------------------
% Files

% create file to store Aggregate_Info.
aggData_fileName = strcat(tracesDirName,'aggSPECTG.mat');
aggData_mFile = matfile(aggData_fileName);

% create file to store CodeWords
cw_fileName = strcat(tracesDirName,'cw_final_SPECTG.mat');
cw_mFile = matfile(cw_fileName);

% create file to store CodeWords-Histograms of each Trace
cwHist_fileName = strcat(tracesDirName,'cwHIST_SPECTG.mat');
cwHist_mFile = matfile(cwHist_fileName,'Writable',true);


% -------------------------------------------------------------------------
% Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
DIMS_SPECTG = aggData_mFile.DIMS_SPECTG ;
CODEWORDS_FINAL = cw_mFile.CODEWORDS;


% -------------------------------------------------------------------------
% Read each song and compute CodeWord Histogram

% --------------------------------------------------------------
    % Determine tau-nearest-neighbor CodeWord_Histogram for each song.
    
    Files = dir(aggregateDataDirName);
    
    countSong = 0;
    for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
        
        if (Files(k).bytes > 0 )
            countSong = countSong + 1
            
            % Read the mat file containing mfcc,dyn_mfcc data
            matfileName=Files(k).name;
            spectg_fileName = strcat(aggregateDataDirName,matfileName);
            spectg_mFile = matfile(spectg_fileName);
            
            
            % extract the data needed.
            DATA = spectg_mFile.NORM_SPECTG;
            [coeffs numFrames] = size(DATA);
            
            
            % find tau-nearest neighbors
            
            size(CODEWORDS_FINAL)
            size(DATA)
            nbrs_of_songFrames = knnsearch(CODEWORDS_FINAL',DATA','k',tau,'distance','euclidean');
            nbrs_of_songFrames = reshape(nbrs_of_songFrames,[],1);
            
            % extract histogram of occurence and re-structure into COL_VEC
            hist_vec = histc(nbrs_of_songFrames,1:numCodeWordsPerBatch);
            hist_vec = reshape(hist_vec,[],1);
            
            % save directly into file
            if countSong > 1
                cwHist_mFile.CODEWORD_HIST(:,countSong) = (1/numFrames) * (1/tau)* hist_vec;
            else
                cwHist_mFile.CODEWORD_HIST = (1/numFrames) * (1/tau)* hist_vec;
            end
            
            clear hist_vec DATA;
            
        end
    end
    
