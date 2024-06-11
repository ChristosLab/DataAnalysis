function outmat = downsample_row(in_mat, n)
%   Quick utility function for downsampling 2D matrices by the rows without LP filtering
inmat = downsample(in_mat', n);
outmat = inmat';
end