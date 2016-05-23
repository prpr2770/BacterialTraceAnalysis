% Script, to generate a field from the data points measured from each cell.
% With an intent to visualize how the fields are evolving with time.
% This visualization technique, might help us understand/visualize any
% particular patterns that's not discernible from the individual
% time-series.

% Script to do the following:
% 1. Iterate over the Ecoli/Styph/Mix files
% 2. From each: generate Fields for each time-step.
% 3. Visualize the motion of the field-membrane.


%%
close all; clear all;

featureType = 'FIELD_VISUALIZATION';%'SPECTOGRAM'; % MFCC

dAnz = 'H:\KraljLab\2016-02-18-PROPS_CALC_Ecoli_vs_Salmonella\ChangeScores';
cd(dAnz)

flist = dir('*.mat');
nfiles = length(flist);

fs = 5; % Sampling frequency. 5Hz

% Directory to save processed data
saveDir = [dAnz filesep featureType '_' num2str(date)];
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
% 1. Read the denoised signals and centroids.
% 2. Visualize the Field Variations.

for fID = 1:nfiles
    fID
    fname = flist(fID).name
    
    % --------------------------------------------------------------------
    % FolderName to be data related to a file.
    fDIR_windows = [saveDir filesep fname(1:end-4)];
    if ~exist(fDIR_windows)
        mkdir(fDIR_windows)
    end
    
    
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
    all_traces = [];
    all_centroids = [];
    
    % iterate over each trace, extract signal and centroids.
    for idx = 1:ncells
        
        % obtain signal
        sig = data(idx).sig_dn_wav;
        centroid = data(idx).centroid;
        
        all_traces = [all_traces; sig];
        all_centroids = [all_centroids; centroid];
        
    end
    % --------------------------------------------------------------------
    % Process the aggregated Information.
    
    % -------------------------------------
    % spatial coordinates estimate.
    min_coords = min(all_centroids);
    max_coords = max(all_centroids);
    
    % ----------------------------------------
    % time-series estimate of the std-deviations/variances.
    
    % compute spatial-mean time-series.
    spat_mean = mean(all_traces,1); % column operation
    
    % compute the deviation of traces
    deviate_all_traces = all_traces - repmat(spat_mean,size(all_traces,1),1);
    
    
    % -----------------------------------------------------------------
    % Time Series Analysis
    
    total_Time = size(all_traces,2);
    
    reference_coords = all_centroids;
    field_values = deviate_all_traces;
    % ------------------------------------------------
    % Visualize the dataset.
    
    fig1 = figure(1);
    colormap parula
    X = reference_coords(:,1);
    Y = reference_coords(:,2);
    Smoothness = 0.00005;
    
    XNodes = linspace(0.5*min(X),1.5*max(X),50);
    YNodes = linspace(0.5*min(Y),1.5*max(Y),50);
    
    minZ = min(min(field_values));
    maxZ = max(max(field_values));
    
    for t=1:total_Time
        t
        Z = field_values(:,t);
        ZNodes = RegularizeData3D(X,Y,Z, XNodes, YNodes, 'interp', 'bicubic', 'smoothness', Smoothness);
        surf(XNodes,YNodes,ZNodes,'facealpha', 0.50)
        zlim([(1.5*minZ) (1.5*maxZ)])
        plt_nm = sprintf('calcField_%d.png',t);
        saveas(fig1, fullfile(fDIR_windows, plt_nm), 'png');
        
        % pause
    end
    
    close all;
end % EndIterateOverFiles





