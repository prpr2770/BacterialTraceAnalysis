% extract approximate motifs, using SAX + Random Projection

close all; clear all;

% comment below line after  1 run.
loadData=0;

if loadData == 0
    try
        tracesDirName = 'H:\KraljLab\isa\';
        tracesFileName = strcat(tracesDirName,'isa_Analysis.mat');
        load(tracesFileName); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

outDirName = 'H:\KraljLab\sax\';

% ------------------------------------------------------------------
% applying sax
real_data = volt_sigs(120,:);
figure;
plot(real_data)


wndw_len = 32;
nseg = 16;
alphabet_size = 4;
[str, ptrs] = timeseries2symbol(real_data, wndw_len, nseg, alphabet_size);

% str - numeric code
% ptr - location of begining of code
plot(str','o')

% % % -------------------------------------------------------------------
% % % visualization:
% % % create multiple strings of same length, which contain these subsequences, concatanated with zeros.
% % % Then do a scatter plot of the differnt coded-subsequences.
% %
% % totalLen = size(str,2)+ size(str,2)-1;
% % subseqs = zeros(size(str,1),totalLen);
% % for idx = 1:size(str,1)
% % subseqs(idx,idx:idx+7) = str(idx,:);
% % end
% %  figure; plot(subseqs')

% -------------------------------------------------------------------
% Random Projection: for a given trace, determine the similarity.
sym_data = str;
projection_dim = nseg/4; % k < nseg - faulty_digits
bucket_threshold = 3;
num_iterations = 32;
%  [adjacency_matrix, similarity_matrix] = getRandProjectionSimilarityMatrix(data, projection_size, bucket_threshold, num_iterations)
getRandProjectionSimilarityMatrix


figure;
subplot(1,2,1)
colormap
imagesc(similarity_matrix)
colorbar
subplot(1,2,2)
colormap
imagesc(adjacency_matrix)
colorbar
% -------------------------------------------------------------------
% Identify the time-domain subsequences that are similar to each other.
% and plot them.

% Inputs: ptrs, real_data, sym_data, adjacency_matrix



for idx = 1:length(adjacency_matrix)
    %      idx = ceil(rand()*length(adjacency_matrix));
    nbrs = find(adjacency_matrix(idx,:));
    nbr_ptrs = ptrs(nbrs);
    
    % after identifying neigbors, collect all the windows corresponding to
    % the nbrs and plot them over the orig data.
    motif_set_mask = zeros(length(nbrs), size(real_data,2));
    motif_subseqs = zeros(length(nbrs), wndw_len);
    for idy = 1:length(nbrs)
        motif_set_mask(idy,nbr_ptrs(idy):nbr_ptrs(idy)+wndw_len -1) = 1;
        motif_subseqs(idy,:) = real_data(nbr_ptrs(idy):nbr_ptrs(idy)+wndw_len -1);
    end
    
    motif_set = motif_set_mask.*repmat(real_data,length(nbrs),1);
    
    % archive it inside structure for all
    all_motifs_struct(idx).motif_set = motif_set;
    all_motifs_struct(idx).motif_subseqs = motif_subseqs;
    
    % plot the motifs over the real-ts
    fig10 = figure(10);
    plot(real_data,'g')
    hold on
    plot(motif_set','r')
    hold off
    fig_title = sprintf('motif subseqs, track: %d',idx);
    title(fig_title)
    
    
    % compare motifs among themselves
    fig20 = figure(20);
    subplot(2,1,1)
    plot(motif_subseqs','Color',[0.8 0.8 0.8]);
    hold on
    mean_motif_subseq = mean(motif_subseqs);
    plot(mean_motif_subseq,'Color',[0.8 0.4 0.8]);
    hold off
    fig_title = sprintf('motifs, track: %d',idx);
    title(fig_title)
    
    % let's normalize the windows and observe them
    subplot(2,1,2)
    subseqs = motif_subseqs';
    subseqs = (subseqs - repmat(mean(subseqs),size(subseqs,1),1))./repmat(std(subseqs),size(subseqs,1),1);
    plot(subseqs,'Color',[0.8 0.8 0.8]);
    hold on
    mean_subseqs = mean(subseqs');
    plot(mean_subseqs,'Color',[0.8 0.4 0.8]);
    hold off
    fig_title = sprintf('z-normalized motifs, track: %d',idx);
    title(fig_title)
end