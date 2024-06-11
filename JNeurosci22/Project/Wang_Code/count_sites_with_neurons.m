%%  No. of sites w/ isolated neurons and w/ delay beta tuning
yl = [0, 60];
counts_tun = zeros(4, 2, 2, 2);
i_range = 1:4;
j_range = 0:1;
k_range = 0:1;
l_range = 0:1;
for i = 1:numel(i_range) %   Monkey
    for j = 1:numel(j_range) %    Stage
        for k = 1:numel(k_range) %   yes_neuron
            for l = 1:numel(l_range) %  Tuning
                counts_tun(i, j, k, l) = sum(...
                    (lfp_sites_tbl.monkey_id == i_range(i)).*...
                    (lfp_sites_tbl.aged == j_range(j)).*...
                    ((lfp_sites_tbl.n_neuron > 0) == k_range(k)).*...
                    (lfp_sites_tbl.delay_beta == l_range(l)).*...
                    (lfp_sites_tbl.task_id == 1)...
                    );
            end
        end
    end
end
figure
to_bar = squeeze(counts_tun(:, 1 ,:, 2)./sum(counts_tun(:, 1, : , :), 4));
bar(to_bar*100)

legend({'w/o neurons', 'w/ neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with delay period beta tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_beta_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/ neurons', 'w/o neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with delay period beta tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_beta_adult'), '-dpng', '-r400');