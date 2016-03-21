% find motifs using moen for the ca_sigs.

close all; clear all;

% comment below line after  1 run.
loadData=0;

if loadData == 0
    try
        tracesDirName = 'H:\KraljLab\isa_noisy\';
        tracesFileName = strcat(tracesDirName,'ca_isa_Analysis.mat');
        load(tracesFileName); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end


outDirName = 'H:\KraljLab\moen\';
% =========================================================================

% idx = ceil(rand()*size(ca_sigs_dn,1));
for idx = 1:size(ca_sigs_dn,1)
    raw_time_series = ca_sigs_dn(idx,:);
    len_TS = length(raw_time_series);
    min_motif_len = 16;
    max_motif_len = len_TS/2;
    
    figure();
    figTrace = gcf;
    plot(raw_time_series);
    % save the fig
    plt_nm = sprintf('motifs_ts_trackID_%d.png',idx);
    saveas(figTrace, fullfile(outDirName, plt_nm), 'png');
    
    % write the raw_time_series into a txt file, to be read by the script.
    FileName = 'ca_raw_ts.txt';
    dlmwrite(FileName,raw_time_series,'\n')
    
    % compute motifs
    % !Moen.exe ca_raw_ts.txt 700 32 400
    c = sprintf('!Moen.exe %s %d %d %d > out.txt',FileName, 700, 32, 128);
    result = evalc(c);
    disp(result);

    %%% At this point, the algorithm is finished, the data has been written to disk
    %%% Ask the user if she wants to see the output in a pretty plot.
    
    L = csvread('out.txt');
    plotMotifs(raw_time_series,L(:,2),L(:,3),L(:,1));
    figMotifs = gcf;
    
    % save the fig
    plt_nm = sprintf('motifs_trackID_%d.png',idx);
    saveas(figMotifs, fullfile(outDirName, plt_nm), 'png');
    
    close all;
end