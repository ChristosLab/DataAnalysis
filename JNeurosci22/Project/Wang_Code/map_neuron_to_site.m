function r = map_neuron_to_site(mapping_mat, sites)
[r, ~] = find(mapping_mat(:, sites));
end