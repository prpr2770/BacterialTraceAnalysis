% Converting Traces into audio files, to use MIR toolbox.
clear all; close all; 
load('H:\KraljLab\PROPS_data.mat')
fs = 5; % samplingFreq 5Hz
trace1 = intens(1,:);

% -------------------------------------------------------------------------

% MIR toolbox
audioFile = miraudio(trace1', 5)                % miraudio(colVec,samplingRate) + samplingRate = 5Hz.
frames = mirframe(audioFile,'Length',20,'s');   % frame-width = 20seconds.
spectra = mirspectrum(frames)
flx = mirflux(spectra)
brite = mirbrightness(frames)                       % time evolution of brightness; higher values => most energy in higher frequencies.
cntr = mircentroid(frames)                          % around which fqs sound energy is centered
simMatrix = mirsimatrix(audioFile,'Frame')


% -------------------------------------------------------------------------
% Auditory-toolbox

% Extract spectrogram
segsize = 128;
nlap = 4;
ntrans = 2;
[specgram,raw] = spectrogram(trace1,segsize,nlap,ntrans);
imagesc(specgram)

% -------------------------------------------------------------------------
% Matlab toolbox: Obtain Spectrogram

sig = trace1';
Nx = length(sig);
nsc = floor(Nx/4.5);
novlap = floor(3*nsc/4);
nfft = max(256,2^nextpow2(nsc));
wind = hamming(nsc);

[SPECTG, F, T] = spectrogram(sig,wind,novlap,nfft,fs)

fig1 = figure(1)
subplot(2,1,1)
%surf(T,F(2:end),abs(SPEC(2:end,:)))
surf(T,F,abs(SPECTG))
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
title 'Regular Plot'

subplot(2,1,2)
row = 2;
surf(T,F(row:end),abs(SPECTG(row:end,:)))
% surf(T,F,abs(SPEC))
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
title 'Regular Plot'

fig2 = figure(2)
imagesc(log10(abs(SPECTG(row:end,:))))


% -------------------------------------------------------------------------
