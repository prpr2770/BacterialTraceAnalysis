% Visualize and compare corresponding motifs
% 1. randomly select idx
% 2. obtain adjacency_matrix col-indices for the voltage-motifs
% 3.

% =========================================================================
close all;
% outDIR = 'H:\KraljLab\studyDenoising\motifs\caDerived\';
inDIR = 'H:\KraljLab\studyDenoising\';
outDIR = [inDIR filesep 'motifs_' num2str(date)];

outDIR = [outDIR filesep 'caDerived'];
if ~exist(outDIR)
    mkdir(outDIR)
end


volt_signatures = [];
ca_signatures = [];

% idx = ceil(rand()*length(ca_motifs));
for idx = 1:length(ca_motifs)
    
    % identify neighbours based on the Ca-Derived ChgScore
    volt_adjacency = ca_motifs_adjacency_matrix_volt(idx,:);
    ca_adjacency = ca_motifs_adjacency_matrix_ca(idx,:);
    
    % Voltage derived motifs
    idx_nbrs_volt = find(volt_adjacency);
    
    nbr_volt_motifs = [];
    
    fig1 = figure(1);
    for id = 1:length(idx_nbrs_volt)
        % for each neighbor
        nbrID = idx_nbrs_volt(id);
        volt_sig = volt_sigs_noisy(nbrID).sig_dn_wav;
        
        % obtain the location of the motifs for Ca/Volt-Derived signals
        leftBorder = ca_motifs(nbrID).leftBorder;
        rightBorder = ca_motifs(nbrID).rightBorder;
        
        motif_mask = zeros(size(volt_sig));
        motif_mask(leftBorder:rightBorder) = 1;
        volt_motif = volt_sig.*motif_mask;
        volt_motif_subseq = volt_sig(leftBorder:rightBorder);
        
        try
            nbr_volt_motifs = [nbr_volt_motifs; volt_motif_subseq];
        catch
            warning('Motif dimension mismatch')
            size(nbr_volt_motifs)
            size(volt_motif)
            
            pause;
            zero_motif = zeros(1,size(nbr_volt_motifs,2));
            nbr_volt_motifs = [nbr_volt_motifs; zero_motif];
            
        end
        % plot/visualize :
        subplot(2,1,1)
        plot(volt_sig,'Color',[0.8,0.8,0.8]);
        hold on;
        subplot(2,1,2)
        plot(volt_motif,'Color',[0.8, 0.4, 0.8]);
        hold on;
        
    end
    
    subplot(2,1,1)
    hold off;
    ylabel('volt')
    subplot(2,1,2)
    hold off;
    ylabel('volt')
    
    set(gcf,'NextPlot','add');
    axes;
    heading = sprintf('Ca-ChgScore derived: similar Voltage motifs for trackID:%d',idx);
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    % -----------------------------------------------------------
    % save the image
    dirName = [outDIR filesep 'voltMotifs'];
    if ~exist(dirName)
        mkdir(dirName)
    end
    plt_nm = sprintf('volt_motifs_trackID_%d.png',idx);
    saveas(fig1, fullfile(dirName, plt_nm), 'png');
    
    % =========================================================================
    % Voltage derived Ca-motifs
    idx_nbrs_ca = find(ca_adjacency);
    
    nbr_ca_motifs = [];
    
    fig2 = figure(2);
    
    for id = 1:length(idx_nbrs_ca)
        % for each neighbor
        nbrID = idx_nbrs_ca(id);
        ca_sig = ca_sigs_noisy(nbrID).sig_dn_wav;
        
        % obtain the location of the motifs for Ca/Volt-Derived signals
        leftBorder = ca_motifs(nbrID).leftBorder;
        rightBorder = ca_motifs(nbrID).rightBorder;
        
        motif_mask = zeros(size(ca_sig));
        motif_mask(leftBorder:rightBorder) = 1;
        ca_motif = ca_sig.*motif_mask;
        ca_motif_subseq = ca_sig(leftBorder:rightBorder);
        % archive
        try
            nbr_ca_motifs = [nbr_ca_motifs; ca_motif_subseq];
        catch
            warning('Motif dimension mismatch')
            size(nbr_ca_motifs)
            size(ca_motif)
            pause;
            
            zero_motif = zeros(1,size(nbr_ca_motifs,2));
            nbr_ca_motifs = [nbr_ca_motifs; zero_motif];
            
        end
        % plot/visualize :
        subplot(2,1,1)
        plot(ca_sig,'Color',[0.8,0.8,0.8]);
        hold on;
        subplot(2,1,2)
        plot(ca_motif,'Color',[0.4, 0.8, 0.8]);
        hold on;
        
    end
    subplot(2,1,1)
    hold off;
    ylabel('calcium')
    subplot(2,1,2)
    hold off;
    ylabel('calcium')
    
    % --------------------------------------------------------------
    set(gcf,'NextPlot','add');
    axes;
    heading = sprintf('Ca-ChgScore derived: similar Calcium motifs for trackID:%d',idx);
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    % save the image
    
    dirName = [outDIR filesep 'caMotifs'];
    %     dirName = strcat(outDIR,'caMotifs\');
    if ~exist(dirName)
        mkdir(dirName)
    end
    plt_nm = sprintf('ca_motifs_trackID_%d.png',idx);
    saveas(fig2, fullfile(dirName, plt_nm), 'png');
    
    % =========================================================================
    % Plot the average volt_motifs and ca_motifs
    
    
    % Process the nbr_volt_motifs and nbr_ca_motifs: Z-normalize each row
    [mean_volt, var_volt, nbr_volt_motifs] = getZNormalize(nbr_volt_motifs);
    [mean_ca, var_ca, nbr_ca_motifs] = getZNormalize(nbr_ca_motifs);
    
    
    
    % obtain the mean-signatures
    nbr_volt_motifs_mean = mean(nbr_volt_motifs,1);
    nbr_ca_motifs_mean = mean(nbr_ca_motifs,1);
    
    try
        volt_signatures = [volt_signatures; nbr_volt_motifs_mean];
        ca_signatures = [ca_signatures; nbr_ca_motifs_mean];
    catch
        warning('Matrix dimensions do not meet.')
        sprintf('volt_signatures  #cols: %d',size(volt_signatures,2))
        sprintf('nbr_volt_motifs_mean  #cols: %d',size(nbr_volt_motifs_mean,2))
        sprintf('ca_signatures  #cols: %d',size(ca_signatures,2))
        sprintf('nbr_ca_motifs_mean  #cols: %d',size(nbr_ca_motifs_mean,2))
        
        pause;
        zero_volt_motif = zeros(1,size(volt_signatures,2));
        zero_ca_motif = zeros(1,size(ca_signatures,2));
        
        volt_signatures = [volt_signatures; zero_volt_motif];
        ca_signatures = [ca_signatures; zero_ca_motif];
        
    end
    
    
    fig3 = figure(3)
    
    subplot(1,2,1) % voltage signature
    plot(nbr_volt_motifs','Color',[0.8 0.8 0.8]);
    hold on
    plot(nbr_volt_motifs_mean,'Color',[0.8 0.4 0.8]);
    hold off
    ylabel('volt')
    
    subplot(1,2,2) % ca signature
    plot(nbr_ca_motifs','Color',[0.8 0.8 0.8]);
    hold on
    plot(nbr_ca_motifs_mean,'Color',[0.8 0.4 0.8]);
    hold off
    ylabel('calcium')
    
    set(gcf,'NextPlot','add');
    axes;
    heading = sprintf('Ca-ChgScore derived: Volt and Ca Motif Signatures for trackID-%d',idx);
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    % save the image
    
    %     dirName = strcat(outDIR,'volt_ca_signatures\');
    dirName = [outDIR filesep 'volt_ca_signatures'];
    if ~exist(dirName)
        mkdir(dirName)
    end
    
    plt_nm = sprintf('signature_motifs_trackID_%d.png',idx);
    saveas(fig3, fullfile(dirName, plt_nm), 'png');
    
    close all;
end

file_nm = 'ca_derived_signatures_motifs_volt_ca.mat';
save(fullfile(outDIR, file_nm),'volt_signatures','ca_signatures');
