function max_class = Neuron_Data_Max(filename)
try
    load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\' filename])
catch
    load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', filename(1:9), '1', filename(10:end)])
end
if ~isempty(MatData)
    if filename(8)=='2'|| filename(8)=='3'
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
    elseif filename(8)=='1'
        for n = 1:8
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
    end
else
    max_class = [];
end
