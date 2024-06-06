function [exam, Classes] = exam_finder_MSNG( filename)

% [max_class_corr, min_class_corr] = MSNG_best_worst(filename);
[max_class_corr, min_class_corr] = MSNG_best_distant(filename);

if rem(max_class_corr,2)==0
    max_classes=[max_class_corr-1, max_class_corr];
elseif rem(max_class_corr,2)~=0 && max_class_corr~=17
    max_classes=[max_class_corr, max_class_corr+1];
elseif max_class_corr==17
    max_classes=max_class_corr;
end
examhigh = ones(1, length(max_classes));

if length( min_class_corr)==1
     min_classes=[min_class_corr, min_class_corr+1];
elseif length(min_class_corr)==2
     min_classes=[min_class_corr(1), min_class_corr(1)+1, min_class_corr(2), min_class_corr(2)+1];
end
examlow = ones(1, length(min_classes)) .* (-1);

exam = [examhigh examlow]; 
Classes = [max_classes min_classes];