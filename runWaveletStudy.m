% Running wavelet study on a single trace instance from the Volt-Ca
% dataset. 

% =======================================================================
% Wavelet Analysis
idx = ceil(rand()*size(volt_sigs,1));

v1 = volt_sigs(idx,:)';
c1 = ca_sigs(idx,:)';

figure;
plot(v1,'r')
hold on
plot(c1,'g')
hold off
legend('volt','calcium','Location','NorthWest')
xlabel('time')
title('Volt-Ca for a given cell')
% -------------------------------------------------------------------
% wavelet denoising
v1_dn = getWaveletDenoisedTrace(v1); 
c1_dn = getWaveletDenoisedTrace(c1);

Nx = length(v1);
figure;
subplot(2,2,1),plot(v1), xlim([1 Nx])
title('Original: volt')
subplot(2,2,3),plot(v1_dn), xlim([1 Nx])
title('De-noised: volt')
subplot(2,2,2),plot(c1), xlim([1 Nx])
title('Original: Ca')
subplot(2,2,4),plot(c1_dn), xlim([1 Nx])
title('De-noised: Ca')


% -------------------------------------------------------------------
% wavelet analysis
fs = 5;

scalesCWT = helperCWTTimeFreqVector(0.01,(fs/2),(fs)/(2*pi),1/fs,32);

cwt_v1 = cwtft({v1,1/5},'wavelet','bump',...
    'scales',scalesCWT,'PadMode','symw');
cwt_c1 = cwtft({c1,1/5},'wavelet','bump',...
    'scales',scalesCWT,'PadMode','symw');

cwt_v1_dn = cwtft({v1_dn,1/5},'wavelet','bump',...
    'scales',scalesCWT,'PadMode','symw');
cwt_c1_dn = cwtft({c1_dn,1/5},'wavelet','bump',...
    'scales',scalesCWT,'PadMode','symw');


tm = 1:length(v1);
figure;
subplot(2,2,1)
helperCWTTimeFreqPlot(cwt_v1.cfs,tm,cwt_v1.frequencies,'surf',...
    'volt','Seconds','Hz')
subplot(2,2,3)
helperCWTTimeFreqPlot(cwt_c1.cfs,tm,cwt_c1.frequencies,'surf',...
    'calcium','Seconds','Hz')

subplot(2,2,2)
helperCWTTimeFreqPlot(cwt_v1_dn.cfs,tm,cwt_v1_dn.frequencies,'surf',...
    'volt-dn','Seconds','Hz')
subplot(2,2,4)
helperCWTTimeFreqPlot(cwt_c1_dn.cfs,tm,cwt_c1_dn.frequencies,'surf',...
    'calcium-dn','Seconds','Hz')


% ----------------------------------------------------------------
% wavelet coherence

% scalesCWT = helperCWTTimeFreqVector(0.01,(fs/2),(fs/2)/(2*pi),1/fs,32);
dt = 1/fs;
scales = helperCWTTimeFreqVector(0.01,(fs/2),fs,dt,32);
wcoh = wcoher(v1,c1,scales./dt,'cmor0.5-5','nsw',8);
wcoh_dn = wcoher(v1_dn,c1_dn,scales./dt,'cmor0.5-5','nsw',8);
frequencies = scal2frq(scales./dt,'cmor0.5-5',dt);

figure;
subplot(2,1,1)
surf(tm,frequencies,abs(wcoh).^2); view(0,90); shading interp;
axis tight;
hc = colorbar;
hc.Label.String = 'Coherence';
title('Wavelet Coherence')
xlabel('Seconds'),ylabel('Hz');

subplot(2,1,2)
surf(tm,frequencies,abs(wcoh_dn).^2); view(0,90); shading interp;
axis tight;
hc = colorbar;
hc.Label.String = 'Coherence';
title('Wavelet Coherence: Denoised')
xlabel('Seconds'),ylabel('Hz');


% ----------------------------------------------------------------
% Error: Test for Gaussianity of the Wavelet Approx Error




