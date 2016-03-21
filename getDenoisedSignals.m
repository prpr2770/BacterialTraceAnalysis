function sigs_dn = getDenoisedSignals(sigs)
% Ensure sigs are row-vectors

sigs = sigs';
sigs_dn = [];
for idx = 1:size(sigs,2)
   sig = sigs(:,idx);
   sig_dn = getWaveletDenoisedTrace(sig); 
   sigs_dn = [sigs_dn sig_dn];
end

% convert into row
sigs_dn = sigs_dn';
end
