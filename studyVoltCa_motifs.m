% Study motif extraction from data.
close all; clear all;

% comment below line after  1 run.
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
    
    %--------------------------------------------------------------------------
    % extract volt and ca signals
    
    totalsigs = size(data,2)
    
    % idx = ceil(rand()*totalsigs);
    volt_sigs = [];
    ca_sigs = [];
    
    for idx = 1:totalsigs
        v1 = data(idx).TimeTraceVolt;
        c1 = data(idx).TimeTraceCa;
        f1 = data(idx).Frame;
        
        if ( idx>1 && (length(v1)~=size(volt_sigs,2) || length(c1)~=size(ca_sigs,2) ))
            dimV = size(volt_sigs,2);
            dimC = size(ca_sigs,2);
            volt_sigs = [volt_sigs; v1(1:dimV)];
            ca_sigs = [ca_sigs; c1(1:dimC)];
            
        elseif (idx>1 && length(v1)==size(volt_sigs,2) && length(c1)==size(ca_sigs,2) )
            volt_sigs = [volt_sigs; v1];
            ca_sigs = [ca_sigs; c1];
            
        elseif idx == 1
            volt_sigs = [volt_sigs; v1];
            ca_sigs = [ca_sigs; c1];
        end
        
        
    end
    volt_sigs = removeMean(volt_sigs);
    ca_sigs = removeMean(ca_sigs);
    
    % study behaviour of mean-volt/ca signals.
    
    mean_volt = mean(volt_sigs,1);
    mean_ca = mean(ca_sigs,1);
    
    figure;
    subplot(2,2,1)
    plot(volt_sigs','Color',[0.8 0.8 0.8])
    hold on
    plot(mean_volt,'Color',[0.4 0.8 0.8])
    hold off
    subplot(2,2,3)
    plot(ca_sigs','Color',[0.8 0.8 0.8])
    hold on
    plot(mean_ca,'Color',[0.4 0.8 0.8])
    hold off
    subplot(2,2,2)
    plot(mean_volt,'Color',[0.4 0.8 0.8])
    subplot(2,2, 4)
    plot(mean_ca,'Color',[0.4 0.8 0.8])
    
end


    % ==================================================================
    % Extract ISA representation for the signals.
    
    % isa_volt = getISA(volt_sigs);
%     isa_ca = getISA(ca_sigs);
%     [isa_sigs, features_win_sigs, eigVectors, eigValues] = getISA(ca_sigs);
%    [KernelMatrix, testSig_win, testSig_features, isa_sigs,  features_win_sigs, codeWords_eigVectors, referenceWords, eigVectors, eigValues] = getISA(ca_sigs);
    getISA
    
    



















% % % % =========================================================================
% % % try
% % %     outDirName = 'H:\KraljLab\getMotifs\';
% % % catch
% % %     warning('Error in creating output directory');
% % % end
% % %

% % % % save/archive in mat file
% % % fn = 'exploreMotifs.mat';
% % % fName = strcat(outDirName,fn);
% % % save(fName,'data');
