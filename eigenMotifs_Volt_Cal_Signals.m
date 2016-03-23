% Column subset selection: For detecting EigenMotifs.

close all; clear all;

inDIR = 'H:\KraljLab\studyDenoising\motifs\voltDerived';

% find the .mat file
cd(inDIR)

flist = dir('*.mat');
nfiles = length(flist);

if nfiles == 1
    % load('H:\KraljLab\studyDenoising\motifs\voltDerived\volt_derived_signatures_motifs_volt_ca.mat');
    fdata = load(flist(1).name);
else
    error('>1 .mat files in directory')
end

% the variables are elements of 1x1 structure fdata.
% extract data
fields = fieldnames(fdata(1));
if numel(fields) ~= 2
    warning('File has more fields than 2.')
end


%% Extract the motifs and visualize

for idx = 1:2
    
    if regexp(fields{idx},'ca_sig')
        sourceType = 'Calcium';
        srcType = 'Ca';
        altType = 'V';
    elseif regexp(fields{idx},'volt_sig')
        sourceType = 'Volt';
        srcType = 'V';
        altType = 'Ca';
    end
    
    
    if idx == 1
        source = fdata(1).(fields{1});
        alt_source = fdata(1).(fields{2});
    else
        source = fdata(1).(fields{2});
        alt_source = fdata(1).(fields{1});
    end
    
    %%
    target = source;
    k = 12;
    [S W]= GeneralizedGreedySelection(target, source, k)
    source_eigMotifs = source(S,:);
    alt_from_source_eigMotifs = alt_source(S,:);
    
    %%
    fig1 = figure(1)
    numPlots = length(S);
    numRows = min(4, length(S));
    numCols = ceil(numPlots/numRows);
    
    for idy =1:length(S)
        subplot(numRows,numCols,idy)
        p1 = plot(source_eigMotifs(idy,:),'r'); hold on;
        p2 = plot(alt_from_source_eigMotifs(idy,:),'g'); hold off;
        if idy == 1
            legend([p1 p2],srcType,altType);
        end
        
    end
    
    
    set(gcf,'NextPlot','add');
    axes;
    heading = ['Volt/Ca Motifs: from ' sourceType ' sigs'];
    h = title(heading);
    set(gca,'Visible','off');
    set(h,'Visible','on');
    
    plt_nm = ['eigMotifs_' fields{idx}];
    saveas(fig1, fullfile(inDIR, plt_nm), 'png');
    
    %%
    close all
end


%{
%% Execute the above code in a more verbose manner.


source = volt_dervied_ca_signatures;
target = source;
k = 10;
[S_ca W]= GeneralizedGreedySelection(target, source, k)
ca_eigMotifs = source(S_ca,:);
volt_from_ca_eigMotifs = volt_dervied_volt_signatures(S_ca,:);

%%

source = volt_dervied_volt_signatures;
target = source;
k = 10;
[S_volt W]= GeneralizedGreedySelection(target, source, k)
volt_eigMotifs = source(S_volt,:);
ca_from_volt_eigMotifs = volt_dervied_ca_signatures(S_volt,:);

%% Visualize Plots:

% volt_from_calcium_motifs
fig1 = figure(1)
numPlots = length(S_ca);
numRows = min(4, length(S_ca));
numCols = ceil(numPlots/numRows);

for idx =1:length(S_ca)
    subplot(numRows,numCols,idx)
    p1 = plot(ca_eigMotifs(idx,:),'r'); hold on;
    p2 = plot(volt_from_ca_eigMotifs(idx,:),'r'); hold off;
    legend([p1 p2],'Ca','V');
    
end


set(gcf,'NextPlot','add');
axes;
heading = ['Volt/Ca Motifs: from Ca-EigMotifs'];
h = title(heading);
set(gca,'Visible','off');
set(h,'Visible','on');

plt_nm = ['ca-eigmotifs'];
saveas(fig1, fullfile(inDIR, plt_nm), 'png');
%}