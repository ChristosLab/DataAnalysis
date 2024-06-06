function max_class = Neuron_Data_Maxcuerate_DistFrom4LOC(Distfilename)
% 2022/11/01 J Zhu
try
    load(Distfilename)
    if length(MatData.class)==4
            try
                for n = 1:4
                    var(n) = mean([MatData.class(n).ntr.cuerate]);
                end
                temp_class = find(var == max(var));
                max_class(1) = temp_class(1);
                max_class(2) = var(temp_class(1));
            catch
                max_class = nan(1,2);                
            end
    else
        disp('wrong dist total classes')
        max_class = nan(1,2);
    end
    
catch
    lasterr
end