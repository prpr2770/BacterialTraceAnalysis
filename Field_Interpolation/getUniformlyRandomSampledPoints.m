function reference_coords = getReferenceFieldPoints(rect_sides,num_points)
% function to generate num_points uniformly sampled from rectangle of given
% side-length. 

reference_coods = rand(num_points,2).*repmat(rect_sides,num_points,1);

end