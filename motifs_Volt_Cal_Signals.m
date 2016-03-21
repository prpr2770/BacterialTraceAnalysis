% Studying motif patterns in all traces, individually ( UNIVARIATE MOTIFS )
% % 1. Obtain location of max-event change.(index)
% % 2. Identify (leftMargin,rightMargin) for the motif from the index.
% % 3. Extract primary motif and save in new structure
% % 4. Z-normalize the primary motif

% Requirements to execute code:
% denoise_Volt_Cal_Signals.m
% visualize_Volt_Cal_Signals.m

% =========================================================================
% 1. Read the original data:

loadData=0;% comment below line after  1 run.

if loadData == 0
    try
        inDIR = 'H:\KraljLab\studyDenoising\';
        fname = strcat(inDIR,'chg_det_volt_ca_sigs.mat');
        load(fname); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end

outDIR = 'H:\KraljLab\studyDenoising\motifs\';
% =========================================================================

% Motif-parameter
leftMargin = 10;
rightMargin = 20;

% SAX
codeWordLength = 8;
alphabetSize = 4;

% SAX - Random Projection
projection_dim = codeWordLength/4; % k < nseg - faulty_digits
num_iterations = 32;
bucket_threshold = 0.6*num_iterations;   % used to determine similarity

% Motif archival
arr_volt_motifs_volt = [];
arr_volt_motifs_ca = [];

arr_ca_motifs_volt = [];
arr_ca_motifs_ca = [];

for idx = 1:length(ca_sigs_noisy)
    % ================================================================
    % Extract signals: both volt and calcium
    % ---------------------------------------------
    volt_sig = volt_sigs_noisy(idx).sig_dn_wav;
    ca_sig = ca_sigs_noisy(idx).sig_dn_wav;
    
    volt_sig_chgScr = volt_sigs_noisy(idx).sig_dn_chg_score;
    ca_sig_chgScr = ca_sigs_noisy(idx).sig_dn_chg_score;
    
    % ----------------------------------------------------------
    % Voltage derived motifs: determine index of max change
    [maxVal, maxIdx] = max(volt_sig_chgScr);
    leftBorder = maxIdx-leftMargin;
    if leftBorder <=0
        leftBorder = 0;
    end
    rightBorder = maxIdx+rightMargin;
    if rightBorder > length(volt_sig)
        rightBorder = length(volt_sig)
    end
    
    % extract motif
    volt_motif = volt_sig(leftBorder:rightBorder);
    volt_motif_sax = timeseries2symbol(volt_motif, length(volt_motif), codeWordLength, alphabetSize);
    
    ca_motif = ca_sig(leftBorder:rightBorder);
    ca_motif_sax = timeseries2symbol(ca_motif, length(ca_motif), codeWordLength, alphabetSize);
    
    % archive motif
    volt_motifs(idx).volt_motif = volt_motif;
    volt_motifs(idx).volt_motif_sax = volt_motif_sax;
    volt_motifs(idx).ca_motif = ca_motif;
    volt_motifs(idx).ca_motif_sax = ca_motif_sax;
    volt_motifs(idx).leftBorder = leftBorder;
    volt_motifs(idx).rightBorder = rightBorder;
    
    arr_volt_motifs_volt = [arr_volt_motifs_volt; volt_motif_sax];
    arr_volt_motifs_ca = [arr_volt_motifs_ca; ca_motif_sax];
    
    % ----------------------------------------------------------
    % Calcium derived motifs: determine index of max change
    [maxVal, maxIdx] = max(ca_sig_chgScr);
    leftBorder = maxIdx-leftMargin;
    if leftBorder <=0
        leftBorder = 0;
    end
    rightBorder = maxIdx+rightMargin;
    if rightBorder > length(ca_sig)
        rightBorder = length(ca_sig)
    end
    
    % extract motif
    volt_motif = volt_sig(leftBorder:rightBorder);
    volt_motif_sax = timeseries2symbol(volt_motif, length(volt_motif), codeWordLength, alphabetSize);
    
    ca_motif = ca_sig(leftBorder:rightBorder);
    ca_motif_sax = timeseries2symbol(ca_motif, length(ca_motif), codeWordLength, alphabetSize);
    
    % archive motif
    ca_motifs(idx).volt_motif = volt_motif;
    ca_motifs(idx).volt_motif_sax = volt_motif_sax;
    ca_motifs(idx).ca_motif = ca_motif;
    ca_motifs(idx).ca_motif_sax = ca_motif_sax;
    ca_motifs(idx).leftBorder = leftBorder;
    ca_motifs(idx).rightBorder = rightBorder;
    
    arr_ca_motifs_volt = [arr_ca_motifs_volt; volt_motif_sax];
    arr_ca_motifs_ca = [arr_ca_motifs_ca; ca_motif_sax];
    
    %{
    % ================================================================
    % Extract signals
    % ---------------------------------------------
    % voltage signals
    sig = volt_sigs_noisy(idx).sig_dn_wav;
    % determine index of max change
    sig_chgScr = volt_sigs_noisy(idx).sig_dn_chg_score;
    [maxVal, maxIdx] = max(sig_chgScr);
    % extract motif
    sig_motif = sig(maxIdx-leftMargin:maxIdx+rightMargin);
    sig_motif_sax = timeseries2symbol(sig_motif, length(sig_motif), codeWordLength, alphabetSize);
    % archive motif
    volt_motifs(idx).motif = sig_motif;
    volt_motifs(idx).motif_sax = sig_motif_sax;
    arr_volt_motifs = [arr_volt_motifs; sig_motif_sax];

        
    % ---------------------------------------------
    % calcium signals
    sig = ca_sigs_noisy(idx).sig_dn_wav;
    % determine index of max change
    sig_chgScr = ca_sigs_noisy(idx).sig_dn_chg_score;
    [maxVal, maxIdx] = max(sig_chgScr);
    % extract motif
    sig_motif = sig(maxIdx-leftMargin:maxIdx+rightMargin);
    sig_motif_sax = timeseries2symbol(sig_motif, length(sig_motif), codeWordLength, alphabetSize);
    % archive motif
    ca_motifs(idx).motif = sig_motif;
    ca_motifs(idx).motif_sax = sig_motif_sax;
    arr_ca_motifs = [arr_ca_motifs; sig_motif_sax];
    
    %}
end

file_nm = 'motifs_volt_ca_sigs.mat';
save(fullfile(outDIR, file_nm),'ca_motifs', 'volt_motifs', 'arr_ca_motifs_volt','arr_ca_motifs_ca', 'arr_volt_motifs_volt','arr_volt_motifs_ca');

% =========================================================================
% Compare similarities between the extracted motifs and Plot Similarity
% Matrix
% -------------------------------------------------------------------------

% Voltage Signals
sym_data = arr_volt_motifs_volt;
getRandProjectionSimilarityMatrix;
volt_motifs_similarity_matrix_volt = similarity_matrix;
volt_motifs_adjacency_matrix_volt = adjacency_matrix;
% -----------------------------
fig1 = figure(1);
subplot(1,2,1)
colormap
imagesc(similarity_matrix)
colorbar
xlabel('similarity matrix')
subplot(1,2,2)
colormap
imagesc(adjacency_matrix)
xlabel('adjacency matrix')
colorbar
% ------------------------------
% set title
set(gcf,'NextPlot','add');
axes;
heading = sprintf('Volt-ChgScore derived Volt-Motifs');
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');
% ------------------------------
plt_nm = sprintf('volt_motifs_adj_matrix_volt.png');
saveas(fig1, fullfile(outDIR, plt_nm), 'png');

% -----------------------------
% -----------------------------
sym_data = arr_volt_motifs_ca;
getRandProjectionSimilarityMatrix;
volt_motifs_similarity_matrix_ca = similarity_matrix;
volt_motifs_adjacency_matrix_ca = adjacency_matrix;


fig2 = figure(2);
subplot(1,2,1)
colormap
imagesc(similarity_matrix)
colorbar
xlabel('similarity matrix')
subplot(1,2,2)
colormap
imagesc(adjacency_matrix)
xlabel('adjacency matrix')
colorbar
% ------------------------------
% set title
set(gcf,'NextPlot','add');
axes;
heading = sprintf('Volt-ChgScore derived Ca-Motifs');
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');
% -----------------------------
plt_nm = sprintf('volt_motifs_adj_matrix_ca.png');
saveas(fig2, fullfile(outDIR, plt_nm), 'png');

close all;
% =========================================================================
% Calcium Signals
sym_data = arr_ca_motifs_volt;
getRandProjectionSimilarityMatrix;
ca_motifs_similarity_matrix_volt = similarity_matrix;
ca_motifs_adjacency_matrix_volt = adjacency_matrix;

fig1 = figure(1);
subplot(1,2,1)
colormap
imagesc(similarity_matrix)
colorbar
xlabel('similarity matrix')
subplot(1,2,2)
colormap
imagesc(adjacency_matrix)
xlabel('adjacency matrix')
colorbar
% ------------------------------
% set title
set(gcf,'NextPlot','add');
axes;
heading = sprintf('Ca-ChgScore derived Volt-Motifs');
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');
% ------------------------------

plt_nm = sprintf('ca_motifs_adj_matrix_volt.png');
saveas(fig1, fullfile(outDIR, plt_nm), 'png');
% -----------------------------
% -----------------------------
sym_data = arr_ca_motifs_ca;
getRandProjectionSimilarityMatrix;
ca_motifs_similarity_matrix_ca = similarity_matrix;
ca_motifs_adjacency_matrix_ca = adjacency_matrix;

fig2 = figure(2);
subplot(1,2,1)
colormap
imagesc(similarity_matrix)
colorbar
xlabel('similarity matrix')
subplot(1,2,2)
colormap
imagesc(adjacency_matrix)
xlabel('adjacency matrix')
colorbar
% ------------------------------
% set title
set(gcf,'NextPlot','add');
axes;
heading = sprintf('Ca-ChgScore derived Ca-Motifs');
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');
% ------------------------------

plt_nm = sprintf('ca_motifs_adj_matrix_ca.png');
saveas(fig2, fullfile(outDIR, plt_nm), 'png');

% =========================================================================
file_nm = 'compare_motifs_volt_ca_sigs.mat';
save(fullfile(outDIR, file_nm),'ca_motifs_similarity_matrix_ca','ca_motifs_adjacency_matrix_ca','ca_motifs_similarity_matrix_volt' ,'ca_motifs_adjacency_matrix_volt','volt_motifs_similarity_matrix_ca','volt_motifs_adjacency_matrix_ca','volt_motifs_similarity_matrix_volt' ,'volt_motifs_adjacency_matrix_volt');

% =========================================================================
visualize_Motifs_Volt_Cal_Signals
