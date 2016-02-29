% function [adjacency_matrix, similarity_matrix] = getRandProjectionSimilarityMatrix(sym_data, projection_dim, bucket_threshold, num_iterations)
% implement random projection method.


total_data_points = size(sym_data,1);

data_dim = size(sym_data,2);
% totalLen = 5; projection_size = 3;
mask = zeros(1,data_dim);
mask(1:projection_dim) = 1;

similarity_matrix = zeros(total_data_points);

% iterate for the total-counts
for iter = 1:num_iterations
    p = randperm(length(mask));
    mask = (mask(p)==1); % create boolean
    
    % create projections
    data_proj = sym_data(:,mask);
    
    % convert projections into binary, concatanate them, and convert into
    % numeric value
    bin_proj = dec2bin(data_proj');
    bin_proj_mat = reshape(bin_proj',[],size(data_proj,1));
    bin_proj_mat = bin_proj_mat';
    proj_numeric = bin2dec(bin_proj_mat);
    
    % compare the unique values
    unique_proj_numerics = unique(proj_numeric);
    
    if length(unique_proj_numerics) ~= length(proj_numeric)
        % hashing-collision has happened.
        for idx = 1:total_data_points
        
            for idy = idx:total_data_points
               if proj_numeric(idx) == proj_numeric(idy)
                   similarity_matrix(idx,idy) =  similarity_matrix(idx,idy) + 1;
                   similarity_matrix(idy,idx) =  similarity_matrix(idy,idx) + 1;
               end
            end
        end
        
    else
        % no collision detected.
        warning('no collision detected')
    end
    
end % iterations

adjacency_matrix = (similarity_matrix >= bucket_threshold);

% end