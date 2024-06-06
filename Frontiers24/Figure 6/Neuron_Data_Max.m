function max_class = Neuron_Data_Max(filename)
load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\', filename])
% disp(filename)
if ~isempty(MatData)
    if filename(8)=='1'
        Classes=[1,2,3,4,5,6,7,8];
    elseif filename(8)=='2' || filename(8)=='3'
        Classes=[1,6];
    end

    for n = Classes
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

