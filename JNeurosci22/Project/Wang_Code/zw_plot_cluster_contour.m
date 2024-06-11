function out_vertices = zw_plot_cluster_contour(clusters_label)
%%ZW_PLOT_CLUSTER_CONTOUR takes a matrix of cluster labels and 
labels = unique(clusters_label(find(clusters_label)));
out_vertices = cell(1, numel(labels));
for i = 1:numel(labels)
    cluster_flag = double(~~(clusters_label == labels(i)));
    vertices = contourc(cluster_flag, [0.5, 0.5]);
    if vertices(2,1) ~= size(vertices, 2) - 1
        fprintf('countourc failed to identify cluster %d', labels(i));
        return
    end
    %   Removes the summary column of the contour matrix
    vertices = vertices(:, 2:end);
    vertices_new = vertices;
    vert_diff = diff(vertices, 1, 2);
    n_old = size(vertices, 2);
    counter_old = 1;
    counter_new = 1;
    while counter_old < n_old
        if any(vert_diff(:, counter_old) == 0)
        else
            if abs(round(vertices(1, counter_old)) - vertices(1, counter_old)) < 0.001
                add_vertice = [vertices(1, counter_old + 1); vertices(2, counter_old)];
            else
                add_vertice = [vertices(1, counter_old); vertices(2, counter_old + 1)];
            end
            vertices_new = [vertices_new(:, 1:counter_new), add_vertice, vertices_new(:, counter_new + 1:end)];
            counter_new = counter_new + 1;
        end
        counter_old = counter_old + 1;
        counter_new = counter_new + 1;
    end
    out_vertices{i} = vertices_new;
end
end