function max_class = Neuron_Data_Maxcuerate_ProFrom8LOC(Profilename)
% 2022/11/01 J Zhu
try
    load(Profilename)
    if length(MatData.class)==8
        if ~isempty(MatData)
            try
                for n = 1:8
                    var(n) = mean([MatData.class(n).ntr.cuerate]);
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