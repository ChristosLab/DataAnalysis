function dif_n = get_relative_rad(tar_class, ref_class)
for i = 1:numel(tar_class)
    for j = 1:numel(ref_class)
        x = tar_class(i) - ref_class(j);
        dif_n(i, j) = x + 8*(x < 0) - 8*((x > 4) || ((x < 0) && (x > -4)));
    end
end
end