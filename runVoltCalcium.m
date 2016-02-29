% Main script calling all future scripts to be run on batch-data
close all;

loadData=0;
if loadData == 0
    try
        tracesDirName = 'H:\KraljLab\';
        tracesFileName = strcat(tracesDirName,'data_voltage_calcium.mat');
        load(tracesFileName); Fs = 5; % samplingFreq 5Hz
        
        loadData=1;
    catch
        warning('Error Reading Input Data');
    end
end


try
    outDirName = 'H:\KraljLab\voltCal_results\';
catch
    warning('Error in creating output directory');
end


%--------------------------------------------------------------------------
% extract volt and ca signals

totalsigs = size(data,2)

% idx = ceil(rand()*totalsigs);

for idx = 1:totalsigs
    v1 = data(idx).TimeTraceVolt;
    c1 = data(idx).TimeTraceCa;
    f1 = data(idx).Frame;
    
    % ----------------------------------------------------------------------
    %run changeScore
    leftMargin = 100;
    
    
    try
        
        fig1 = figure(1);
        [v1_score_fin, v1_pks, v1_loc] = getChangeScore(v1,leftMargin);
    catch
        warning('Error in extracting Change-Score: VOLT')
    end
    
    try
        fig2 = figure(2);
        [c1_score_fin, c1_pks, c1_loc] = getChangeScore(c1,leftMargin);
    catch
        warning('Error in extracting Change-Score: CALCIUM')
    end
    
    fig3 = figure(3)
    subplot(2,1,1)
    plot(v1_score_fin,'r')
    subplot(2,1,2)
    plot(c1_score_fin,'g')
    
    try
        % savefig as image
        fn = sprintf('volt_trackID_%d.png',idx);
        plotFileName = strcat(outDirName,fn);
        saveas(fig1, fullfile(outDirName, fn),'png');
        
        fn = sprintf('calc_trackID_%d.png',idx);
        plotFileName = strcat(outDirName,fn);
        saveas(fig2, fullfile(outDirName, fn),'png');
        
        fn = sprintf('chgScore_trackID_%d.png',idx);
        plotFileName = strcat(outDirName,fn);
        saveas(fig3, fullfile(outDirName, fn),'png');
    catch
        warning('Error in SAVING PLOTS.')
    end
    
    
    data(idx).voltChangeScore = v1_score_fin;
    data(idx).calcChangeScore = c1_score_fin;
    
    data(idx).calcPeakLocs = c1_loc;
    data(idx).voltPeakLocs = v1_loc;
    close all;
end



fn = 'changeScore_Analysis.mat';
fName = strcat(outDirName,fn);
save(fName,'data');
        
% 
% 
% data(idx).voltChangeScore = v1_score_fin;
% data(idx).calcChangeScore = c1_score_fin;
% 
% data(idx).calcPeakLocs = c1_loc;
% data(idx).voltPeakLocs = v1_loc;
% v1 = data(idx).TimeTraceVolt;
% c1 = data(idx).TimeTraceCa;
% f1 = data(idx).Frame;
