function getInstFreq_OfTracePair(Traces, delays, x1_delayed, x2_delayed,Fs)
% compute Instantaneous Freq and detect synchronization.

% ------------------------------------------------------------------
% zero-delay signal
z1 = hilbert(Traces(:,1));
z2 = hilbert(Traces(:,2));

% ratio of InstFreqs
f1 = Fs/(2*pi)*diff(unwrap(angle(z1)));
f2 = Fs/(2*pi)*diff(unwrap(angle(z2)));
synchRatio = f1./f2;
InstFreqs = [f1 f2];

% relative phase/ generalized phase difference
phi1 = unwrap(angle(z1));
phi2 = unwrap(angle(z2));

%deterministic sense
phi_12 = phi1 - phi2;
figure;
subplot(1,2,1)
plot(phi_12)
subplot(1,2,2)
plot(phi2,phi1)
set(gcf,'NextPlot','add');
axes;
h = title('Relative Phase: Deterministic');
set(gca,'Visible','off');
set(h,'Visible','on');



%statistical sense
psi_12 = mod(phi_12,2*pi);
figure;
h = histogram(psi_12);
title('Relative Phase: Statistical')

% plots
figure;
subplot(1,3,1)
plot(Traces)
title('Traces')

subplot(1,3,2)
plot([f1 f2])
title('InstFreq')

subplot(1,3,3)
plot(synchRatio)
title('Ratio of InstFreq')

% ------------------------------------------------------------------
% delay signals: Relative Phase
Z1_delayed = hilbert(x1_delayed);
Z2_delayed = hilbert(x2_delayed); 

phi1_delayed = unwrap(angle(Z1_delayed));
phi2_delayed = unwrap(angle(Z2_delayed));

f1_delayed = Fs/(2*pi)*diff(phi1_delayed);
f2_delayed = Fs/(2*pi)*diff(phi2_delayed);

% Ratio of instFreq
synchRatio_delayedF1_F2 = f1_delayed./repmat(f2,1,length(delays));
synchRatio_F1_delayedF2 = repmat(f1,1,length(delays))./f2_delayed;

% Relative Phase
relPhase_delayedF1_F2 = phi1_delayed - repmat(phi2,1,length(delays));
relPhase_F1_delayedF2 = repmat(phi1,1,length(delays))- phi2_delayed;



% plots : SYNCH-RATIO
figure;
subplot(1,2,1)
plot(synchRatio_delayedF1_F2)
subplot(1,2,2)
plot(synchRatio_F1_delayedF2)
set(gcf,'NextPlot','add');
axes;
h = title('Ratio of InstFreqs: Delayed');
set(gca,'Visible','off');
set(h,'Visible','on');


figure;
histogram(synchRatio_F1_delayedF2)
title('SynchRatio: delayedX2 + X1')

figure;
histogram(synchRatio_delayedF1_F2)
title('SynchRatio: delayedX1 + X2')

% plots : REL-PHASE
figure;
histogram(mod(relPhase_delayedF1_F2,2*pi))
title('Rel-Phase: delayedX1 + X1')

figure;
histogram(mod(relPhase_F1_delayedF2,2*pi))
title('Rel-Phase: delayedX2 + X1')

figure;
subplot(1,2,1)
plot(relPhase_delayedF1_F2)
subplot(1,2,2)
plot(relPhase_F1_delayedF2)
set(gcf,'NextPlot','add');
axes;
h = title('RELATIVE PHASE: Delayed');
set(gca,'Visible','off');
set(h,'Visible','on');


% ------------------------------------------------------------------



end