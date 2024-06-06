function [exam, Classes] = exam_finder_2_classes( filename,r1r2)

max_class_corr = Neuron_Data_Max(filename);
if filename(8)=='2'
    if max_class_corr(1) == 1
        Classes=[1 2 3 4 5  6  7  8  9 10]; % for remember 1st
        exam   =[1 1 1 1 1 -1 -1 -1 -1 -1];
    elseif max_class_corr(1) ==6
        Classes=[6 7 8 9 10 1  2  3  4  5]; %for remember 1st
        exam =  [1 1 1 1 1 -1 -1 -1 -1 -1];
    else
        Classes = [];
        exam=[];
    end


elseif filename(8)=='1'
    if max_class_corr(1) <=4
        Classes = [max_class_corr(1) max_class_corr(1)+4]; %for remember 1st
        exam=[1 -1];
    elseif max_class_corr(1) >4 && max_class_corr(1) <=8
        Classes = [max_class_corr(1) max_class_corr(1)-4];%for remember 1st
        exam=[1 -1];
    else
        Classes = [];
        exam=[];
    end
else 
    Classes = [];
        exam=[];
end
