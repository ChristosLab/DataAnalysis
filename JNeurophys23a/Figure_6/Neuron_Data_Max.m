function max_class = Neuron_Data_Max(filename)
load(filename)
if ~isempty(MatData)
    for n = [1,6]
       try
           if ~isempty(MatData.class(n).ntr)
               var(n) = mean([MatData.class(n).ntr.cuerate]);
           end
       catch
           disp('Class does not exist')
       end

    end
    temp_class = find(var == max(var));
    max_class(1) = temp_class(1);
    max_class(2) = var(temp_class(1));    
else
    max_class = [];
end
