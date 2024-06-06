function [max_class, min_class] = MSNG_best_worst(filename)
load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\' filename])
if ~isempty(MatData)
    for n = 1:length(MatData.class)
       try
           if ~isempty(MatData.class(n).ntr)
               var(n) = mean([MatData.class(n).ntr.cuerate]);
           end
       catch
           disp('Class does not exist')
       end

    end
    temp_class_max = find(var == max(var));
    max_class(1) = temp_class_max(1);
    max_class(2) = var(temp_class_max(1));  

    temp_class_min = find(var == min(var));
    min_class(1) = temp_class_min(1);
    min_class(2) = var(temp_class_min(1));  
    
else
    max_class = [];
    min_class = [];
end
