% Extract ISA representation for the signals.
% function [KernelMatrix, testSig_win, testSig_features, isa_sigs,  features_win_sigs, codeWords_eigVectors, referenceWords, eigVectors, eigValues] = getISA(sigs)
% sig: row-vectors of time-series.
% % Input:
% % sigs : row vectors of time-series
% %
% % Output:
% % isa_sigs                : time-series sigs rept as isa_features
% % features_win_sigs       : windowed components of each signal, rept by features.
% % codeWords_eigVectors    : final vectors to describe projections.
% % referenceWords          : words used to extract kernel.
% % eigVectors              : all eigvectors of the algorithm
% % eigValues               : all eigvalues of algorithm
% % % ===================================================================


% UNIT-TEST - criteria:
sigs = ca_sigs; % row-data-vectors.

% Parameters
W = 32; % num-samples per window ( sigs are sampled at 5Hz )
delay = W/4;
K = 16; %k-nearest nbrs for graph-laplacian
kai = 1;
Dprime = W;

maxN = 8192; % total number of code-words to be considered.

tracesDirName = 'H:\KraljLab\isa\';

% windows of i'th cell: win_sigs(t,cellID,windowID);
[win_sigs, winCounts] = getWindows(sigs, W, delay);

% extract freq-domain features for the windows of the cells.
features_win_sigs = getFeatureRepresentation(win_sigs,W);     % bag-of-words. - col-vectors

clear win_sigs;

save('isa_Analysis.mat','features_win_sigs','win_sigs','ca_sigs','volt_sigs')

% obtain graph-laplacian - numClusters(K)

[Lnorm, SimilarityMatrix, KernelMatrix, referenceWords] = getGraphLaplacian(features_win_sigs, K, maxN);
save('isa_graph_Analysis.mat','Lnorm', 'KernelMatrix', 'SimilarityMatrix','referenceWords')

% =========================================================================
% solving eigenvalue problem
A = eye(size(KernelMatrix)) + kai*Lnorm*KernelMatrix;
B = KernelMatrix;
[eigVectors,eigValues] = eig(A,B);
size(eigVectors)

[Y,IDS] = sort(diag(eigValues),'ascend');
size(Y)
size(IDS)

% extract subset of all eigVectors we are interested in.
codeWords_eigVectors = eigVectors(:,IDS(1:Dprime));
figure;
colormap hot
imagesc(codeWords_eigVectors)
colorbar
title('code-word-eigVectors')


% ----------------------------------------------------------------
% =========================================================================
% ----------------------------------------------------------------
% Analyze for all songs

close all;
N = size(referenceWords,2);


for testSig = 1:size(sigs,1)
    
    % obtain bag-of-feature-words.
    
    testSig_features = features_win_sigs(:,testSig,:);
    testSig_features = reshape(testSig_features,size(testSig_features,1),[]);
    
    % --------------------------------------
    % create into 3D to subtract from all avail-window frames
    M = size(testSig_features,2);
    
    tSig_win_mat = repmat(testSig_features,1,1,N);
    
    refW_mat = reshape(referenceWords,size(referenceWords,1),1,[]);
refW_mat = repmat(refW_mat,1,M,1);
size(refW_mat)
    % --------------------------------------
    
    diff_mat = (tSig_win_mat - refW_mat);
    dist_mat = sum(diff_mat.^2,1);% 2D-matrix
    
    
    % compute Kernel
    sigma = 2;
    kernel_win = exp(-dist_mat.^2/(2*sigma^2));
    kernel_win = reshape(kernel_win, size(kernel_win,2), [],1) ; % [MxN]
    
    %multiply with Dprime words
    testSig_isa = kernel_win*codeWords_eigVectors;
    
    fig1 = figure(1);
    colormap hot   
    subplot(1,2,1)
    imagesc(testSig_features)
    colorbar
    title('short-time fourier')
    subplot(1,2,2)
    caxis([-1,0])
    imagesc(testSig_isa')
    colorbar
    title('isa')
    % test
    isa_sigs = rand(1,2);
    
    
    fn = sprintf('ca_trackID_%d.png',testSig);
    plotFileName = strcat(tracesDirName,fn);
    saveas(fig1, fullfile(tracesDirName, fn), 'png');
    close all
    
    
end


% end