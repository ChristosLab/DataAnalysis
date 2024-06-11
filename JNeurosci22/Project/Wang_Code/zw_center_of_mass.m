function out = zw_center_of_mass(in_array)
%1-D center of mass
pos = 1:numel(in_array);
out = sum(in_array.*pos)/sum(in_array);
end