function imf_mat = getIMF_matrix(imf)
% imf: structure contining row-vectors of imf

imf_mat = [];
numIMFs = length(imf);
for i=1:numIMFs
    imf_mat = [imf_mat; imf{i}] ;
end

imf_mat = imf_mat';
end