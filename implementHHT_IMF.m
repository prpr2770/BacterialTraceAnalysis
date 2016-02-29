% Implementing huang-Hilbert Transform
% -------------------------------------------------------------------------

close all;clear all;


tracesDirName = 'H:\KraljLab\';
tracesFileName = strcat(tracesDirName,'Data_20-Jan-2016.mat');
load(tracesFileName);
Fs = 5; % samplingFreq 5Hz

ORG_DATA = data;
% extract the values
tt=[]; for i = 1:size(ORG_DATA,2) tt = [tt; ORG_DATA(i).TimeTrace]; end
cntrs=[]; for i = 1:size(ORG_DATA,2) cntrs = [cntrs; ORG_DATA(i).Centroid]; end

% datafile
A = tt;

% % % extract imfs for each trace and save into data structure
% % for i = 1:size(ORG_DATA,2) 
% %     trace = ORG_DATA(i).TimeTrace;
% %     imf = emd(trace);
% %     
% %     % obtain col-vec form of imf
% %     imf_mat = [];
% %     numIMFs = length(imf);
% %     for i=1:numIMFs
% %         imf_mat = [imf_mat; imf{i}] ;
% %     end
% %     imf_mat = imf_mat';
% % Z = hilbert(imf_mat);% % 
% % instFreq = Fs/(2*pi)*diff(unwrap(angle(Z)));

% %     % store into structure
% %     ORG_DATA(i).HILBERT_IMFs = Z;
% %     ORG_DATA(i).INST_FQs = instFreq;
% %     % compute the instantaneous frequencies


% % end


% -------------------------------------------------------------------------
trace = tt(3,:);
imf = emd(trace);

% combine imfs from structure to matrix
imf_mat = [];
numIMFs = length(imf);
for i=1:numIMFs
    imf_mat = [imf_mat; imf{i}] ;
end
imf_mat = imf_mat';

figure; plot(imf_mat)

%-----------------------------------------------------
% compare the reconstruction
reconstruct = sum(imf_mat(:,3:numIMFs)');
err_residue = trace - reconstruct;

figure;
subplot(2,1,1)
plot(reconstruct,'g')
hold on
plot(trace,'r')
hold off
subplot(2,1,2)
plot(err_residue)
%-----------------------------------------------------

Z = hilbert(imf_mat);
instFreq = Fs/(2*pi)*diff(unwrap(angle(Z)));
figure;
plot(instFreq)


t_seq = (1:length(trace))*(1/Fs);
figure;
plot(t_seq(2:end),instFreq)
xlabel('Time')
ylabel('Hz')
grid on
title('Instantaneous Frequency')

figure;
numPlots = size(instfreq,2);
for k=1:3
subplot(3,1,k)
if k~=3
plot(t_seq(2:end),instfreq(:,k))
else
plot(t_seq(2:end),sum(instfreq(:,3:end)'))
end
grid on
end

% -------------------------------------------------------------------------
% determine PHASE-CORRELATIONS between bacteria



%--------------------------------------------------------------------------
% detect Events in a trace
sig = tt(566,:);

B = sig;
alpha = .0;
n = 50;
k = 10;

score1 = change_detection(B,n,k,alpha);
score2 = change_detection(B(:,end:-1:1),n,k,alpha);

subplot(2,1,1);
plot(B, 'b-', 'linewidth',1);
axis([-inf,size(B,2),-inf,inf])
title('Original Signal')

subplot(2,1,2);
score2 = score2(end:-1:1); % why do this?

% 2*n+k-2 is the size of the "buffer zone".
plot([zeros(1,2*n-2+k),score1 + score2], 'r-', 'linewidth',1);
axis([-inf,size(B,2),-inf,inf])
title('Change-Point Score')

