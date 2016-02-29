function norm_sig = getNormalized(sig)

% sig = sig - min(sig);
% norm_sig = sig/sum(sig);

sig = sig - mean(sig);
norm_sig = sig/norm(sig,2);
end