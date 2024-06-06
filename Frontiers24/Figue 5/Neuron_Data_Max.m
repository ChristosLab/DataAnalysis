function max_class = Neuron_Data_Max(filename)
try
    load(['ALLDataCorrErr\' filename])
catch
    load(['ALLDataCorrErr\', filename(1:9), '1', filename(10:end)])
end
if ~isempty(MatData)
    if filename(8)=='2'|| filename(8)=='3'
        for n = [1,6]
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
        temp_class = find(var == max(var));
        max_class(1) = temp_class(1);
        max_class(2) = var(temp_class(1));
        [~, max_class(3)] = baseline_firingrate(MatData);
        [~, max_class(4)] = cuedelay_firingrate(MatData, max_class(1));
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
        [~, max_class(3)] = baseline_firingrate(MatData);
        [~, max_class(4)] = cuedelay_firingrate(MatData, max_class(1));
    end
else
    max_class = [];
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

