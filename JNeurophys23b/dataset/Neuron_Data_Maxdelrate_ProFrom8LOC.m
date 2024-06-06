function max_class = Neuron_Data_Maxdelrate_ProFrom8LOC(Profilename)
% J Zhu, 2022-03-12, for ODR(ver. 2019)
try
    load(Profilename)
    if ~isempty(MatData)
        if length(MatData.class)==8
            try
                for n = 1:8
                    var(n) = mean([MatData.class(n).ntr.cuedelay]);
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