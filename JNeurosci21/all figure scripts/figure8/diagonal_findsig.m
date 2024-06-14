%get data for origin
load('post_decodematch_orimean.mat')
load('post_decodematch_ori95.mat')
load('post_decodematch_all.mat')
load('post_decodematch_mean.mat')
x=linspace(-1,5,57);
xq=linspace(-1,5,400);
interp_95=smooth(interp1(x,confidence_95,xq),10);

%get timepoint for sig difference
x=linspace(-1,5,57);
xq=linspace(-1,5,400);
win_size=20;
win_step=10;
load('pre_113_all.mat');
%dataA=all_results;
dataA=store_mat;
for i=1:size(dataA,1)
    temp_data=dataA(i,:);
    new_A(i,:)=smooth(interp1(x,temp_data,xq),10);
end
load('pre_213_all.mat');
%dataB=all_results;
dataB=store_mat;
for i=1:size(dataB,1)
    temp_data=dataB(i,:);
    new_B(i,:)=smooth(interp1(x,temp_data,xq),10);
end
sig_time=find_sigtime(new_A,new_B,win_size,win_step);
%find start and end point
if length(sig_time)>1
    find_discontinue=diff(sig_time);
    start_index=[1];
    end_index=length(sig_time);
    if length(find(find_discontinue>1))>0
       start_index=[start_index,find(find_discontinue>1)+1];   
       end_index=[find(find_discontinue>1),end_index]; 
    end
%calculate start and end sigtime in real world time
real_start=-1+(sig_time(start_index)-1)*6/400;
real_end=-1+(sig_time(end_index)-1)*6/400;
end
