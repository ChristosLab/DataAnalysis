function max_class = Neuron_Data_Maxcuerate_ProFrom8LOC(datain)
% 2022/11/06 J Zhu
try
    load(Profilename)
    if length(MatData.class)==8
        if ~isempty(MatData)
            try
                for n = 1:8
                    if isfield(MatData.class(n).ntr,'fixrate') % diff data file format
                        var(n) = mean([MatData.class(n).ntr.cuerate] - [MatData.class(n).ntr.fixrate]);
                    else
                        var(n) = mean([MatData.class(n).ntr.cuerate] - [MatData.class(n).ntr.fix]);
                    end
                end
                temp_class = find(var == max(var));
                max_class(1) = temp_class(1);
                max_class(2) = var(temp_class(1));
            catch
                max_class = nan(1,2);
            end
        else
            max_class = nan(1,2);
        end
    else
        disp('wrong prosac total classes')
    end

catch
    lasterr
end