% Supervised Clustering: PROPS_ [Ecoli + Styph]
% Compute: CodeWordHistogram_ITML_TSNE


%%
close all; clear all;

dAnz = 'H:\KraljLab\2016-03-17 PROPS in E coli vs salmonella\NoisyResults_18-Mar-2016\SetA\ChangeScores_20-Mar-2016\spectG_21-Mar-2016';
cd(dAnz)

fs = 5; % Sampling frequency. 5Hz

% create file to store CodeWords-Histograms of each Trace
cwHist_fileName = [dAnz filesep 'cwHIST_SPECTG.mat'];
cwHist_mFile = matfile(cwHist_fileName,'Writable',true);

%% Iterate through the sub-folders

dirList = dir();

folderCount = 0;
for dirID = 1:length(dirList)
    dirName = dirList(dirID).name
    if (length(regexp(dirName,'styph')) || length(regexp(dirName,'ecoli')))
        %%
        
        folderCount = folderCount + 1;
        
        % Navigate into folder
        dirLoc = [dAnz filesep dirName];
        cd(dirLoc)
        flist = dir('NORM_SPECTG_*.mat');
        nfiles = length(flist);
        
        
        % Compile all NORM_SPECTG
        ALL_NORM_SPECTG = []; % Col-Vectors
        for fid = 1:nfiles
            fdata = load(flist(fid).name);
            ALL_NORM_SPECTG = [ALL_NORM_SPECTG fdata(1).NORM_SPECTG];
        end
        
        % Archive the NORM_SPECTG
        FILE_NORM_SPECTG(folderCount).ALL_NORM_SPECTG = ALL_NORM_SPECTG;
        FILE_NORM_SPECTG(folderCount).FID_NORM_SPECTG = folderCount*ones(1,size(ALL_NORM_SPECTG,2));
        
    else
        
        warning('No ECOLI and SALMONELLA DataFolders Available in DIR> ')
        
    end
    
end


%% Compile all the NORM_SPECTG into one array.

GLOBAL_NORM_SPECTG = [];
GLOBAL_NORM_SPECTG_clusterID = [];

for fid = 1:length(FILE_NORM_SPECTG)
    GLOBAL_NORM_SPECTG = [GLOBAL_NORM_SPECTG FILE_NORM_SPECTG(fid).ALL_NORM_SPECTG];
    GLOBAL_NORM_SPECTG_clusterID = [GLOBAL_NORM_SPECTG_clusterID FILE_NORM_SPECTG(fid).FID_NORM_SPECTG];
end

%% Use K-Means clustering to detect Code-Words.

numCodeWords = 256;
iter = 20;                      % for knn-ma-toolbox.

CODEWORDS = ma_kmeans(GLOBAL_NORM_SPECTG', iter, numCodeWords);
CODEWORDS = CODEWORDS'; % Convert into Col-vectors.

cwHist_mFile.CODEWORDS = CODEWORDS;

%% Compute CODEWORD_histogram for each trace of each folder

tau = 5;
GLOBAL_CW_HIST = [];
GLOBAL_CW_HIST_clusterID = [];

folderCount = 0;
for dirID = 1:length(dirList)
    dirName = dirList(dirID).name
    if (length(regexp(dirName,'styph')) || length(regexp(dirName,'ecoli')))
        %%
        
        folderCount = folderCount + 1;
        
        % Navigate into folder
        dirLoc = [dAnz filesep dirName];
        cd(dirLoc)
        flist = dir('NORM_SPECTG_*.mat');
        nfiles = length(flist);
        
        %% iterate over the files
        for fid = 1:nfiles
            fdata = load(flist(fid).name);
            DATA = fdata(1).NORM_SPECTG;     % col-vectors
            numFrames = size(DATA,2);
            
            nbrs_of_Frames = knnsearch(CODEWORDS',DATA','k',tau,'distance','euclidean');
            nbrs_of_Frames = reshape(nbrs_of_Frames,[],1); % col-vector 1D
            
            % extract histogram of occurence and re-structure into COL_VEC
            hist_vec = histc(nbrs_of_Frames,1:numCodeWords);
            hist_vec = reshape(hist_vec,[],1);
            hist_vec = (1/numFrames) * (1/tau)* hist_vec;
            
            % archive the code-word histogram for each trace
            GLOBAL_CW_HIST = [GLOBAL_CW_HIST hist_vec];
            GLOBAL_CW_HIST_clusterID = [GLOBAL_CW_HIST_clusterID folderCount];
        end
        
    else
        
        warning('No ECOLI and SALMONELLA DataFolders Available in DIR> ')
        
    end
    
end


cwHist_mFile.GLOBAL_CW_HIST = GLOBAL_CW_HIST;
cwHist_mFile.GLOBAL_CW_HIST_clusterID = GLOBAL_CW_HIST_clusterID;


%% Unsupervised Clustering: tSNE - Visualization


% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 30;

figure;
tsne_DATA = tsne(GLOBAL_CW_HIST', GLOBAL_CW_HIST_clusterID', out_dims, initial_dims, perplexity);
PLOT_TITLE = sprintf('tSNE Embedding of CodeWordsHistogram-Spectogram');

% -------------------------------------------------------------------------
% scatter plot for all the data
fig = figure()
colormap jet
h = scatter3(tsne_DATA(:,1),tsne_DATA(:,2),tsne_DATA(:,3),'filled');
title(PLOT_TITLE)

figName = [dAnz filesep 'cwHist_tsne_data.fig'];
savefig(figName);


%% Supervised Clustering: ITML + tSNE-Visualization
[dist_metric, dist_matrix] = runITMLonDataSet(GLOBAL_CW_HIST', GLOBAL_CW_HIST_clusterID');

% implement tSNE on Distance Matrix
perplexity = 30;
out_dims = 3;

figure;
itml_tsne_DATA = tsne(dist_matrix, GLOBAL_CW_HIST_clusterID', out_dims, perplexity);

clusterKeys = unique(GLOBAL_CW_HIST_clusterID);
% scatter plot for all the data
fig13 = figure(13)
colormap(parula(length(clusterKeys)));

for i = 1:length(clusterKeys)
% for i = 2:length(clusterKeys)-1
    cellsOfACluster = find(GLOBAL_CW_HIST_clusterID == i);
    hue = (length(clusterKeys)+1-i)*ones(length(cellsOfACluster),1);
    rad = 30*ones(length(cellsOfACluster),1);
    h = scatter3(itml_tsne_DATA(cellsOfACluster,1),itml_tsne_DATA(cellsOfACluster,2),itml_tsne_DATA(cellsOfACluster,3),rad,hue,'filled')
    hold on
end
legend(num2str(clusterKeys))
hold off

figName = [dAnz filesep 'cwHist_itml_tsne_data.fig'];
savefig(figName);
