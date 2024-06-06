function [max_class, max_class_cuerate, fixrate, max_class_cdrate, distantclass] = MSNG_best_distant(filename)
load(['C:\Users\mozumdr\OneDrive - Vanderbilt\Desktop\Backup\Lab\SophiaProject_FS_RS\ALLDataCorrErr\' filename])
if ~isempty(MatData)
    for n = 1:length(MatData.class)
       try
           if ~isempty(MatData.class(n).ntr)
               var(n) = mean([MatData.class(n).ntr.cuerate]);
           else
               var(n) = nan;
           end
       catch
           var(n) = nan;
           disp('Class does not exist')
       end

    end
    temp_class_max = find(var == max(var));
    max_class = temp_class_max(1); 
    max_class_cuerate = var(max_class);
    [~, fixrate] = baseline_firingrate(MatData);
    [~, max_class_cdrate]  = cuedelay_firingrate(MatData, max_class);
%     if isnan(max_class_cuerate) || isnan(fixrate) || isnan(max_class_cdrate)
%         pause
%     end

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
end

function [all_fixrates, fixrate]    = baseline_firingrate(MatData)

    fixrate      = [];
    all_fixrates = [];


    if ~isempty(MatData)
        for c = 1:length(MatData.class)
            if ~isempty(MatData.class(c).ntr)
                all_fixrates = [all_fixrates, mean([MatData.class(c).ntr.fix])];
            else
                disp('Empty class!!')
            end
        end
        fixrate = mean(all_fixrates);
    else
        disp('Empty Matdata')
    end
end

function [all_delayrates, cuedelayrate]    = cuedelay_firingrate(MatData, max_class)

    cuedelayrate      = [];
    all_delayrates = [];


    if ~isempty(MatData)
        for c = max_class
            if ~isempty(MatData.class(c).ntr)
                all_delayrates = [all_delayrates, [MatData.class(c).ntr.cuedelay]];
            else
                disp('Empty class!!')
            end
        end
        cuedelayrate = mean(all_delayrates);
    else
        disp('Empty Matdata')
    end
end
