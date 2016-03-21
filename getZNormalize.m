function     [mean_sig, var_sig, z_sigs] = getZNormalize(sigs)
% z-Normalize the row-vectors
    mean_sig = mean(sigs,2);
    mean_sig = repmat(mean_sig,1,size(sigs,2));
    var_sig = sum((sigs - mean_sig).^2,2);
    var_sig = repmat(var_sig,1,size(sigs,2));
    z_sigs = (sigs - mean_sig )./var_sig;
    
end