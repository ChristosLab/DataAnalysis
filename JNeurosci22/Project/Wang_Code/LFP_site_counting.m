sum(lfp_sites_tbl.Pro_presac_250ms_.*lfp_sites_tbl.ODR_delay_gamma_tuning)/sum(lfp_sites_tbl.Pro_presac_250ms_)
sum(ismember(lfp_sites_tbl.task_id, 'ODR').*(lfp_sites_tbl.Pro_presac_250ms_ == 0).*lfp_sites_tbl.ODR_delay_gamma_tuning)/sum(ismember(lfp_sites_tbl.task_id, 'ODR').*(lfp_sites_tbl.Pro_presac_250ms_ == 0))
%%
sum(lfp_sites_tbl.Pro_C.*lfp_sites_tbl.ODR_delay_gamma_tuning)/sum(lfp_sites_tbl.Pro_presac_250ms_)
sum(ismember(lfp_sites_tbl.task_id, 'ODR').*(lfp_sites_tbl.Pro_C == 0).*lfp_sites_tbl.ODR_delay_gamma_tuning)/sum(ismember(lfp_sites_tbl.task_id, 'ODR').*(lfp_sites_tbl.Pro_presac_250ms_ == 0))
%% No. of sites with neurons; No. of sites with responsive neurons
counts_res = zeros(4, 2, 2);
i_range = 1:4;
j_range = 0:1;
k_range = 0:1;
for i = 1:numel(i_range) %   Monkey
    for j = 1:numel(j_range) %    Stage
        for k = 1:numel(k_range) %   Responsive
            counts_res(i, j, k) = sum((lfp_sites_tbl.monkey_id == i_range(i)).*(lfp_sites_tbl.aged == j_range(j)).*(lfp_sites_tbl.responsive == k_range(k)).*(lfp_sites_tbl.task_id == 1));
        end
    end
end
%%  No. of sites with neurons; No. of sites with isolated neurons
counts_neuron = zeros(4, 2);
i_range = 1:4;
j_range = 0:1;
for i = 1:numel(i_range) %   Monkey
    for j = 1:numel(j_range) %    Stage
        counts_neuron(i, j) = sum((lfp_sites_tbl.monkey_id == i_range(i)).*(lfp_sites_tbl.aged == j_range(j)).*(lfp_sites_tbl.n_neuron > 0).*(lfp_sites_tbl.task_id == 1));
    end
end
%%
for i = 1:4
    i
    disp('res')
    counts_res(i,:,2)
    disp('neuron')
    counts_neuron(i, :)
end
%%  No. of sites w/ responsive neurons and w/ delay gamma tuning
counts_tun = zeros(4, 2, 2, 2);
i_range = 1:4;
j_range = 0:1;
k_range = 0:1;
l_range = 0:1;
for i = 1:numel(i_range) %   Monkey
    for j = 1:numel(j_range) %    Stage
        for k = 1:numel(k_range) %   Responsive
            for l = 1:numel(l_range) %  Tuning
                counts_tun(i, j, k, l) = sum(...
                    (lfp_sites_tbl.monkey_id == i_range(i)).*...
                    (lfp_sites_tbl.aged == j_range(j)).*...
                    (lfp_sites_tbl.responsive == k_range(k)).*...
                    (lfp_sites_tbl.delay_gamma == l_range(l)).*...
                    (lfp_sites_tbl.n_neuron > 1).* ...
                    (lfp_sites_tbl.task_id == 1)...
                    );
            end
        end
    end
end
figure
to_bar = squeeze(counts_tun(:, 1 ,:, 2)./sum(counts_tun(:, 1, : , :), 4));
bar(to_bar)
set(gca, 'XTickLabel', monkey_id)
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar)
set(gca, 'XTickLabel', monkey_id)
%%  No. of sites w/ isolated neurons and w/ delay gamma tuning
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
                    (lfp_sites_tbl.delay_gamma == l_range(l)).*...
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
title('Percent of LFP sites with delay period gamma tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_gamma_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/ neurons', 'w/o neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with delay period gamma tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_gamma_adult'), '-dpng', '-r400');
%%  No. of sites w/ isolated neurons and w/ cue gamma tuning
yl = [0,60]
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
                    (lfp_sites_tbl.cue_gamma == l_range(l)).*...
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
title('Percent of LFP sites with cue period gamma tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_gamma_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/o neurons', 'w/ neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with cue period gamma tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_gamma_adult'), '-dpng', '-r400');
%%  No. of sites w/ cue responsive neurons and w/ cue gamma tuning (Not good)
yl = [0,60]
counts_tun = zeros(4, 2, 2, 2);
i_range = 1:4;
j_range = 0:1;
k_range = 0:1;
l_range = 0:1;
for i = 1:numel(i_range) %   Monkey
    for j = 1:numel(j_range) %    Stage
        for k = 1:numel(k_range) %   Pro_Cue
            for l = 1:numel(l_range) %  Tuning
                counts_tun(i, j, k, l) = sum(...
                    (lfp_sites_tbl.monkey_id == i_range(i)).*...
                    (lfp_sites_tbl.aged == j_range(j)).*...
                    (lfp_sites_tbl.Pro_C == k_range(k)).*...
                    (lfp_sites_tbl.cue_gamma == l_range(l)).*...
                    (lfp_sites_tbl.n_neuron > 1).* ...
                    (lfp_sites_tbl.task_id == 1)...
                    );
            end
        end
    end
end
figure
to_bar = squeeze(counts_tun(:, 1 ,:, 2)./sum(counts_tun(:, 1, : , :), 4));
bar(to_bar*100)

legend({'w/o cue response', 'w/ cue response'}, 'Location', 'northwest')
title('Percent of LFP sites with cue period gamma tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'cue_resp_cue_gamma_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/o cue response', 'w/ cue response'}, 'Location', 'northwest')
title('Percent of LFP sites with cue period gamma tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'cue_resp_cue_gamma_adult'), '-dpng', '-r400');
%%  No. of sites w/ isolated neurons and w/ cue beta tuning
yl = [0,60]
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
                    (lfp_sites_tbl.cue_beta == l_range(l)).*...
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
title('Percent of LFP sites with cue period beta tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_beta_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/o neurons', 'w/ neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with cue period beta tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_beta_adult'), '-dpng', '-r400');
%%  No. of sites w/ isolated neurons and w/ cue alphabeta tuning
yl = [0,60]
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
                    (lfp_sites_tbl.cue_alphabeta == l_range(l)).*...
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
title('Percent of LFP sites with cue period alphabeta tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_alphabeta_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/o neurons', 'w/ neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with cue period alphabeta tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_cue_alphabeta_adult'), '-dpng', '-r400');
%%  No. of sites w/ isolated neurons and w/ delay alphabeta tuning
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
                    (lfp_sites_tbl.delay_alphabeta == l_range(l)).*...
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
title('Percent of LFP sites with delay period alphabeta tuning (adolescent)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_alphabeta_young'), '-dpng', '-r400');
% 
figure
to_bar = squeeze(counts_tun(:, 2 ,:, 2)./sum(counts_tun(:, 2, : , :), 4));
bar(to_bar*100)

legend({'w/ neurons', 'w/o neurons'}, 'Location', 'northwest')
title('Percent of LFP sites with delay period alphabeta tuning (adult)')
set(gca, 'XTickLabel', monkey_id)
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1);
ylim(yl)
ytickformat('percentage')
print(gcf, fullfile(project_dir, fig_lib, 'neuron_delay_alphabeta_adult'), '-dpng', '-r400');