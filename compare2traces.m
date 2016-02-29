% What can you tell about two traces?

function compare2traces(traces,Fs)
% function compare2traces
% ---------------------------------------
% Initialize data to test function
% % % traceLen = 900;
% % % traces = rand(traceLen,2);
% % % Fs = 5;

figure;
plot(traces)
title('Original Traces')

% emd - analysis
x1 = traces(:,1); 
x2 = traces(:,2);
imf_x1 = emd(x1);
imf_x2 = emd(x2);
imf_x1 = getIMF_matrix(imf_x1);
imf_x2 = getIMF_matrix(imf_x2);

% % % determine best representative IMF: least difference norm
% % x1_numIMFs = size(imf_x1,2);
% % x2_numIMFs = size(imf_x2,2);
% % 
% % residue_1 = x1 - imf_x1(:,end);
% % residue_2 = x2 - imf_x2(:,end);
% % 
% % diff_x1 = imf_x1 - repmat(residue_1,1,x1_numIMFs);
% % diff_x2 = imf_x2 - repmat(residue_2,1,x2_numIMFs);
% % 
% % err_x1 = (sum(diff_x1.*diff_x1)).^(0.5);
% % err_x2 = (sum(diff_x2.*diff_x2)).^(0.5);
% % 
% % [val, x1_indx] = min(err_x1)
% % [val, x2_indx] = min(err_x2)

% % % using best rept-IMF: averaging middle IMFs
% % x1_indx = 5; x2_indx = 5;
% % x1_est = sum(imf_x1(:,3:5)')';
% % x2_est = sum(imf_x1(:,3:5)')';

x1_est = imf_x1(:,end);
x2_est = imf_x2(:,end);
traces = [x1_est x2_est];

% get vel-acc: visualize
[normVel, normAcc] = getVelAcc_OfTracePair(traces,Fs);

% obtain autoCorrelation and Cross-Corelation functions
[Traces, delays, x1_delayed, x2_delayed, corr_x1_dx1, corr_x1_dx2, corr_x2_dx1, corr_x2_dx2] = getCorr_OfTracePair(traces,Fs);

% obtain instantaneous phases
 getInstFreq_OfTracePair(Traces, delays, x1_delayed, x2_delayed,Fs);

end
