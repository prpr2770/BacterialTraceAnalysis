function PhaseSynchIndexMatrix = getPhasSynchronizationMatrix(tt,Fs)

traces = tt'; %column vector

estTraces = [];
for idx = 1:size(traces,2)
imf = emd(traces(:,idx));
imf = getIMF_matrix(imf);
est = sum(imf(:,3:end-2)')'; % removing first IMF component as noise
estTraces = [estTraces est];
end

traces = estTraces(10:end-10,:);

figure;
subplot(1,2,1)
plot(tt)
title('Orig Traces')
subplot(1,2,2)
plot(traces)
title('EMD Traces')

Z = hilbert(traces);
PHI = unwrap(angle(Z));

n=1; m=1;
numTraces = size(traces,2);
PhaseSynchIndexMatrix = [];
for idx = 1:numTraces
    phi = PHI(:,idx);
    relPhase = n*repmat(phi,1,numTraces)- m*PHI;
    INDX = (mean(cos(relPhase)).^2 + mean(sin(relPhase)).^2 ).^(0.5);    
    PhaseSynchIndexMatrix = [PhaseSynchIndexMatrix; INDX];
end


threshold = 0.90;
PhaseSynchIndexMatrix = getClusterSortedMatrix(PhaseSynchIndexMatrix, threshold);
figure;
imagesc(PhaseSynchIndexMatrix);
figure;
imagesc(PhaseSynchIndexMatrix>threshold);




end