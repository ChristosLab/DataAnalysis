function center_coord = zw_find_cluster_center(labelmat)
%FIND_CLUSTER_CENTER computes the center according to the 4-direction
%limits a connected region
center_coord = round([mean(find(any(labelmat, 2))), mean(find(any(labelmat, 1)))]);
end