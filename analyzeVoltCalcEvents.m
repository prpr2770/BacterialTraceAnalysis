% analyse volt-calcium event waveforms. to extract information about the 
% voltage induced calcium events. 


close all;
% clear all;

loadData = 1;
if loadData
    try
        inDirName = 'H:\KraljLab\voltCal_results\';
        fName = strcat(inDirName ,'changeScore_Analysis.mat');
        load(fName);
        totalsigs = size(data,2);
        loadData = 0;
    catch
        error('Data loading error.')
    end
end


% -------------------------------------------------------------------------
%  identifying common events 
% given a refernce time-index, consider the interval (-20, +40)
leftBorder = 20;
rightBorder = 40;

% grpA_intervals = {};
% grpB_intervals = {};

idx = ceil(rand()*totalsigs);
countA = 0;
countB = 0;
idx = 0;
% for idx = 1:2%totalsigs
while(idx < totalsigs)
    
    idx = idx + 1;
    
    V = data(idx).TimeTraceVolt;
    C = data(idx).TimeTraceCa;
    
    changeV = data(idx).voltChangeScore;
    changeC = data(idx).calcChangeScore;
    
    % ------------------------------------------------------------------------
    % Group A
    t_volt_groupA = data(idx).voltProlongedEventLocs;
    t_calc_groupA = data(idx).calcProlongedEventLocs;
    
    % extract volt-cal events from trace, only if events exist.
    if length(t_volt_groupA) || length(t_calc_groupA)
        [volt_intervals, volt_subseq_TS] = getSubsequence(V,t_volt_groupA,leftBorder,rightBorder);
        [calc_intervals, calc_subseq_TS] = getSubsequence(C,t_calc_groupA,leftBorder,rightBorder);
        countA = countA + 1;
        
        grpA_intervals{countA}.voltEvents =  volt_intervals;
        grpA_intervals{countA}.calcEvents =  calc_intervals;
        grpA_intervals{countA}.traceID = idx;
        grpA_intervals{countA}.Volt = V;
        grpA_intervals{countA}.Calc = C;
        grpA_intervals{countA}.volt_subseq_TS = volt_subseq_TS;
        grpA_intervals{countA}.calc_subseq_TS = calc_subseq_TS;
        
        % save plots f1 and f2
        %         delete fig1 fig2;
        
        % Plot all subsequences in this
        fig5 = figure(5);
        ax1 = subplot(2,1,1);
        % plot voltage intervals
        N = size(volt_intervals,2);
        for k =1:N
            sig = V;
            plot(sig ,'Color',[0.80,0.80,0.80]);
            hold on
            sig = volt_intervals{1,k}.sig;
            time = volt_intervals{1,k}.time;
            plot(time,sig ,'Color',[0.80,0.0,0.40]);
            hold on
        end
        ylabel('volt')
       ylim(ax1,[0 3000]);
        %         subplot(2,2,3)
        %         plot(changeV,'Color',[0.80,0.80,0.80]);
        %         hold on
        %
        
        ax2= subplot(2,1,2);
        % plot calcium intervals
        N = size(calc_intervals,2);
        for k =1:N
            sig = C;
            plot(sig ,'Color',[0.80,0.80,0.80]);
            ylim([0 6000])
            hold on
            sig = calc_intervals{1,k}.sig;
            time = calc_intervals{1,k}.time;
            plot( time, sig, 'Color',[0.80,0.0,0.40]);
            hold on
        end
        ylabel('calcium')
        ylim(ax2,[0 6000]);
        title('volt-triggered events')

        %         subplot(2,2,4)
        %         plot(changeC,'Color',[0.80,0.80,0.80]);
        %         hold on
        
    end
    
    % ------------------------------------------------------------------------
    % ========================================================================
    % ------------------------------------------------------------------------
    
    % ------------------------------------------------------------------------
    % Group B
    t_volt_groupB = data(idx).voltImpulseEventLocs;
    t_calc_groupB = data(idx).calcImpulseEventLocs;
    
    % extract volt-cal events from trace, only if events exist.
    if length(t_volt_groupB) || length(t_calc_groupB)
        [volt_intervals, volt_subseq_TS] = getSubsequence(V,t_volt_groupB,leftBorder,rightBorder);
        [calc_intervals, calc_subseq_TS] = getSubsequence(C,t_calc_groupB,leftBorder,rightBorder);
        countB = countB + 1;
        
        grpB_intervals{countB}.voltEvents =  volt_intervals;
        grpB_intervals{countB}.calcEvents =  calc_intervals;
        grpB_intervals{countB}.traceID = idx;
        grpB_intervals{countB}.Volt = V;
        grpB_intervals{countB}.Calc = C;
        grpB_intervals{countB}.volt_subseq_TS = volt_subseq_TS;
        grpB_intervals{countB}.calc_subseq_TS = calc_subseq_TS;
        
        % save plots f1 and f2
        %         delete fig1 fig2;
        
        % Plot all subsequences in this
        fig10 = figure(10);
        ax1 = subplot(2,1,1);
        % plot voltage intervals
        N = size(volt_intervals,2);
        for k =1:N
            sig = V;
            plot(sig ,'Color',[0.80,0.80,0.80]);
            hold on
            sig = volt_intervals{1,k}.sig;
            time = volt_intervals{1,k}.time;
            plot(time,sig ,'Color',[0.80,0.0,0.40]);
            hold on
        end
        ylabel('volt')
        ylim(ax1,[0 3000]);
        %         subplot(2,2,3)
        %         plot(changeV,'Color',[0.80,0.80,0.80]);
        %         hold on
        %
        
        ax2 = subplot(2,1,2);
        % plot calcium intervals
        N = size(calc_intervals,2);
        for k =1:N
            sig = C;
            plot(sig ,'Color',[0.80,0.80,0.80]);
            hold on
            sig = calc_intervals{1,k}.sig;
            time = calc_intervals{1,k}.time;
            plot( time, sig, 'Color',[0.80,0.0,0.40]);
            hold on
        end
        ylabel('calcium')
        ylim(ax2,[0 6000]);
        title('NON-volt-triggered events')
    end
    
    % ------------------------------------------------------------------------a
end


fName = strcat(inDirName ,'volt-calcium-events.mat');
    
save(fName,'grpB_intervals','grpA_intervals');