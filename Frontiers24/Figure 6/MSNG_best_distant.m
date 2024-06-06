function [max_class, distantclass] = MSNG_best_distant(filename)
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
    max_class = temp_class_max(1); 

    if max_class<9
        distantclass = 15;
    elseif max_class>8 && max_class<17
        distantclass = 7;
    elseif max_class == 17
        distantclass = [15, 7];
    end
    
else
    max_class = [];
    distantclass = [];
end
