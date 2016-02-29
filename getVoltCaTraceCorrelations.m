% extract information regarding the correlation between the calcium and
% voltage events. 

% Tasks:
% 1. plot histogram of the intervals between the volt-calcium events. 
% 2. Plot histogram of max time at which Ca or V events occured. 
% 3. Autocorrelation between the signals. 

clear all;
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

outDirName = 'H:\KraljLab\voltCal_results\xcorr\interCell\';

tDelays_all = zeros(1,totalsigs);
tDelays_chg = zeros(1,totalsigs);

xcorr_VC_sum = [];
xcorr_VC_CHG_sum = [];
count = 0;

for idx = 1:totalsigs
    idx
chgV = data(idx).voltChangeScore;
chgC = data(idx).calcChangeScore;
V = data(idx).TimeTraceVolt;
C = data(idx).TimeTraceCa;

[~,maxV_indx] = max(chgV);
[~,maxC_indx] = max(chgC);


% % % % +==============================================================
% extract MEAN and then Normalize. 
V = getNormalized(V);
C = getNormalized(C);

chgV = getNormalized(chgV);
chgC = getNormalized(chgC);

[xcorr_VC,lag] = xcorr(V,C);

[maxVal, tDelay] = max(xcorr_VC);
tDelays_all(idx) = lag(tDelay);

[xcorr_VC_CHG,lag_CHG] = xcorr(chgV,chgC);
[maxVal, tDelay] = max(xcorr_VC_CHG);
tDelays_chg(idx) = lag_CHG(tDelay);

if(idx>1 && (length(xcorr_VC_sum) == length(xcorr_VC)) &&  (length(xcorr_VC_CHG) == length(xcorr_VC_CHG_sum)))
xcorr_VC_sum = [xcorr_VC_sum;  xcorr_VC];
xcorr_VC_CHG_sum = [xcorr_VC_CHG_sum; xcorr_VC_CHG];
count = count + 1;
elseif idx ==1
xcorr_VC_sum = [xcorr_VC_sum;  xcorr_VC];
xcorr_VC_CHG_sum = [xcorr_VC_CHG_sum; xcorr_VC_CHG];
count = count + 1;
lag_VC_sum = lag;
lag_VC_CHG_sum = lag_CHG;
end

% % % ---------------------------- 
% % fig1 = figure(1)
% % subplot(2,2,1)
% % plot(V,'r')
% % hold on
% % plot(C,'g')
% % hold off
% % ylabel('volt(R), Calcium(G)')
% % 
% % subplot(2,2,3)
% % plot(chgV,'r')
% % hold on
% % plot(chgC,'g')
% % hold off
% % ylabel('chgV(R), chgC(G)')
% % 
% % 
% % subplot(2,2,2)
% % plot(lag,xcorr_VC)
% % title('xcorr V-C')
% % 
% % subplot(2,2,4)
% % plot(lag_CHG,xcorr_VC_CHG)
% % title('xcorr chgV-chgC')
% % % ------------------------------
% % % save plot figure
% % 
% % fn = sprintf('trackID_%d.png',idx);
% % plotFileName = strcat(outDirName,fn);
% % saveas(fig1, fullfile(outDirName, fn), 'png');
% % close all

end

xcorr_VC_mean = sum(xcorr_VC_sum)/count;
xcorr_VC_CHG_mean = sum(xcorr_VC_CHG_sum)/count;

fig1 = figure(1)
subplot(2,1,1)
plot(lag_VC_sum, xcorr_VC_mean)
xlabel('volt-ca xcorr')
subplot(2,1,2)
plot(lag_VC_CHG_sum, xcorr_VC_CHG_mean)
xlabel('volt-ca chg xcorr')
fn = 'avg_xcorr.png';
saveas(fig1, fullfile(outDirName, fn), 'png');


fig2 = figure(2);
subplot(2,1,1)
plot(tDelays_all)
ylabel('V-C')
subplot(2,1,2)
plot(tDelays_chg,'r')
ylabel('chgV-chgC')
title('\tau-delay from cross-correlation')

fn = 'tau-delay.png';
saveas(fig2, fullfile(outDirName, fn), 'png');

% plot the correlations for all signals as an image
fig3 = figure(3);
subplot(2,1,1)
imagesc(xcorr_VC_sum)
colorbar
title('xcorr V-C')
subplot(2,1,2)
imagesc(xcorr_VC_CHG_sum)
colorbar
title('xcorr chgV-chgC')
fn = 'all_xcorr.png';
saveas(fig3, fullfile(outDirName, fn), 'png');

% % % [val, indx] = min(tDelays_all);
% % % 
% % % 
% % % idx = indx;
% % % V = data(idx).TimeTraceVolt;
% % % C = data(idx).TimeTraceCa;
% % % 
% % % C2 = circshift(C,[0 100]);
% % % [crossCorr,lag] = xcorr(C2,V);
% % % [maxVal, tDelay] = max(crossCorr);
% % % 
% % % % ---------------------------- 
% % % figure;
% % % subplot(3,1,1)
% % % plot(V)
% % % ylabel('volt')
% % % subplot(3,1,2)
% % % plot(C)
% % % hold on
% % % plot(C2)
% % % hold off
% % % ylabel('Ca')
% % % subplot(3,1,3)
% % % plot(lag,crossCorr)
% % % ylabel('cross correlation')
% % % 
