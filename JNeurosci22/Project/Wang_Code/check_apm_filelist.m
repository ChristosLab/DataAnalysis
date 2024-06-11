f_list = [prosac_files; antisac_files];
fname_list = cell(numel(prosac_files) + numel(antisac_files), 1);
for i = 1:numel(fname_list)
    if regexp(f_list(i).name, '.*catch.*')
        fname_list{i} = '';
        continue
    end
    fname_list{i} = f_list(i).name(1:8);
end
%%
fname_list = unique(fname_list);
fname_list = fname_list(2:end);
out_tbl = cell2table(fname_list, 'VariableNames', {'FileName'});
%%
writetable(out_tbl, fullfile(project_dir, output_database, 'third_batch_files.csv'));
%%
d1 = dir('Z:\Mack\files for Zhengyang\LEM');
d2 = dir('Z:\Mack\files for Zhengyang\IND');
d3 = dir('Z:\Mack\files for Zhengyang\KEN');
d4 = dir('Z:\Mack\files for Zhengyang\JAC');
%%
ds = [d1; d2; d3; d4];
%%
added_name ={ds.name};
for i = 1:numel(added_name)
    try
    added_name{i} = lower(added_name{i}(1:8));
    catch
        continue
    end
end