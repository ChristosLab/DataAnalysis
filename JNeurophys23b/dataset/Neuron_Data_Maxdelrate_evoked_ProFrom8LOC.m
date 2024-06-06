function max_class = Neuron_Data_Maxdelrate_evoked_ProFrom8LOC(Profilename)
% J Zhu, 20221106, for ODR(ver. 2019)
try
    load(Profilename)
    if ~isempty(MatData)
        if length(MatData.class)==8
            try
                for n = 1:8
                    if isfield(MatData.class(n).ntr,'fixrate') % diff data file format
                        var(n) = mean([MatData.class(n).ntr.cuedelay] - [MatData.class(n).ntr.fixrate]);
                    else
                        var(n) = mean([MatData.class(n).ntr.cuedelay]- [MatData.class(n).ntr.fix]);
                    end
                end
                temp_class = find(var == max(var(1:8)));
                max_class(1) = temp_class(1);
                max_class(2) = var(temp_class(1));
            catch
                max_class = nan(1,2);
            end
        else
            disp('wrong ODR total classes')
            max_class = nan(1,2);
        end
    else
        disp('Empty MatData')
        max_class = nan(1,2);
    end
catch
    lasterr
    max_class = nan(1,2);
end