%find significant cells 
[num,txt,raw]=xlsread('neuron_location.xlsx',3);
neuron_num=cell2num([raw(2:end,3)]);
neuron_loc=string([raw(2:end,4)]);
PFC_area={'46','8'};
Parietal_are={'LIP','PPC','7a'};
pfc_index=[];
parietal_index=[];
location_list=[15,16,13,14,11,12,9,10,17,1,2,3,4,5,6,7,8];
sample_list=[1,5,2,5,3,5,4,5,5,6,5,7,5,8,5,9,5];
for i=1:length(neuron_num)
    if ismember(neuron_loc{i},PFC_area)
    pfc_index=[pfc_index,i];
    else
    parietal_index=[parietal_index,i];
    end
end
pfc_cellnum=neuron_num(pfc_index);
parietal_cellnum=neuron_num(parietal_index);
load('new_msngcorrect_info.mat');
load('new_msngcorrect_data.mat');
all_session=cell_info{1,1};
for i=1:length(all_session)
    temp_sessionname=all_session{i};
    temp_animalname=temp_sessionname(1:3);
    if strcmpi(temp_animalname,'ken')
        all_animalcode(i)=1;
    else
        all_animalcode(i)=2;
    end
end
%animal specific selection

%[C,select_index,ib]=intersect(cell_info{3},pfc_cellnum);%
[C,select_index,ib]=intersect(cell_info{3},parietal_cellnum);
%select_index=intersect(select_index,find(all_animalcode==2));
Parietal_proportion_resample=[0.3065,0.2607];%Parietal main factor proportion
PFC_proportion_resample=[0.4329,0.4192];%PFC main factor proportion
resample_size=200;


good_data=ones(1,length(select_index));
for i=1:length(select_index)

    cell_data=msng_data(select_index(i),:);  
    all_cuedelay=[];
    all_cue=[];
    all_sample=[];
    all_delay=[];
    all_cueloc=[];
    cuedelay_loc=[];
    cuedelay_epoch=[];
    for j=1:length(cell_data)
        cueclass_data=cell_data{location_list(j)}; %use only match trials
        cue_class_trials(j)=length(cueclass_data);
        cueclass_fix_rate=[];
        cueclass_cue_rate=[];
        cueclass_cuedelay_rate=[];
        cueclass_sample_rate=[];
        if cue_class_trials(j)>0
           cueclass_cueon=[cueclass_data.Cue_onT];
           cueclass_sampleon=[cueclass_data.Sample_onT];
           cueclass_TS={cueclass_data.TS};
           for p=1:cue_class_trials(j)
               cueclass_fix_rate(p)=length(find(cueclass_TS{p}>0 & cueclass_TS{p}<cueclass_cueon(p)));
               cueclass_cue_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p) & cueclass_TS{p}<cueclass_cueon(p)+0.5))/0.5;
               cueclass_cuedelay_rate(p)=length(find(cueclass_TS{p}>cueclass_cueon(p)+0.5 & cueclass_TS{p}<cueclass_sampleon(p)))/3;
               cueclass_sample_rate(p)=length(find(cueclass_TS{p}>cueclass_sampleon(p) & cueclass_TS{p}<cueclass_sampleon(p)+0.5))/0.5;            
           end
            
   %     if length(cueclass_ismatch)==length(cueclass_fix_count)

        cue_loc=j*ones(1,length(cueclass_cue_rate));
        sample_loc=sample_list(j)*ones(1,length(cueclass_cue_rate));   
        cuedelay_loc=[cuedelay_loc,cue_loc,cue_loc];
        all_cueloc=[all_cueloc,cue_loc];
        cuedelay_epoch=[cuedelay_epoch,[ones(1,length(cueclass_cue_rate)),2*ones(1,length(cueclass_cue_rate))]];
        all_cuedelay=[all_cuedelay,cueclass_cue_rate,cueclass_cuedelay_rate];
        all_cue=[all_cue,cueclass_cue_rate];
        all_sample=[all_sample,cueclass_sample_rate];
        all_delay=[all_delay,cueclass_cuedelay_rate];
        temp_cuerate(j)=mean(cueclass_cue_rate);
        temp_samplerate(j)=mean(cueclass_sample_rate);
        temp_delayrate(j)=mean(cueclass_cuedelay_rate);
%        end
        end
    end
    %location_to_fit=[-2:1:2];
    location_to_fit=[-2,-2,-1,-1,-0.5,-0.5,-0.25,-0.25,0,0.25,0.25,0.5,0.5,1,1,2,2];

    try
    gauss_cue(i,:)=coeffvalues(fit(location_to_fit',(temp_cuerate-min(temp_cuerate))','gauss1'));
    [obj_cue,gof_cue,opt_cue]=fit(location_to_fit',(temp_cuerate-min(temp_cuerate))','gauss1');
    gof_cue_r2(i)=gof_cue.rsquare;
    gauss_sample(i,:)=coeffvalues(fit(location_to_fit',(temp_samplerate-min(temp_samplerate))','gauss1'));
    [obj_sample,gof_sample,opt_sample]=fit(location_to_fit',(temp_samplerate-min(temp_samplerate))','gauss1');
    gof_sample_r2(i)=gof_sample.rsquare;
    gauss_delay(i,:)=coeffvalues(fit(location_to_fit',(temp_delayrate-min(temp_delayrate))','gauss1'));
    [obj_delay,gof_delay,opt_delay]=fit(location_to_fit',(temp_delayrate-min(temp_delayrate))','gauss1');
    gof_delay_r2(i)=gof_delay.rsquare;
  
    if gauss_cue(i,2)>2
         gauss_cue(i,2)=2;
      elseif gauss_cue(i,2)<-2
         gauss_cue(i,2)=-2;
      end
      if gauss_sample(i,2)>2
         gauss_sample(i,2)=2;
      elseif gauss_sample(i,2)<-2
         gauss_sample(i,2)=-2;
      end
      if gauss_delay(i,2)>2
         gauss_delay(i,2)=2;
      elseif gauss_delay(i,2)<-2
         gauss_delay(i,2)=-2;
      end
    catch
    gauss_cue(i,:)=[0,0,0]; 
    gauss_sample(i,:)=[0,0,0];
    gauss_delay(i,:)=[0,0,0];
    end
      
    if min(cue_class_trials)<5 %filter number of trials
        good_data(i)=0; 
     %   stim_p(i)=[];
     %   delay_p(i)=[];
    end
    cell_numtrials(i)=sum(cue_class_trials);
    cell_maxtrial(i)=max(cue_class_trials);
    
    if length(all_cuedelay)>0 & cue_class_trials(9)>0
    store_cue(i,:)=temp_cuerate;
    store_sample(i,:)=temp_samplerate;
    store_delay(i,:)=temp_delayrate;
    [temp,raw_cue_preferloc(i)]=max(temp_cuerate);
    [temp,raw_sample_preferloc(i)]=max(temp_samplerate);
    [temp,raw_delay_preferloc(i)]=max(temp_delayrate);
    [cuevsdelay_p(i,:),tbl_cuevsdelay]=anovan(all_cuedelay,{cuedelay_loc cuedelay_epoch},'model','interaction','display','off');
    f_cuevsdelay(i,:)=cell2mat(tbl_cuevsdelay(2:4,6))';
    [cue_p(i,:),tbl_cue]=anovan(all_cue,{all_cueloc},'display','off');
    f_cue(i,:)=cell2mat(tbl_cue(2,6));
    [sample_p(i,:),tbl_sample]=anovan(all_sample,{all_cueloc},'display','off');
    f_sample(i,:)=cell2mat(tbl_sample(2,6));
    [delay_p(i,:),tbl_delay]=anovan(all_delay,{all_cueloc},'display','off');
    f_delay(i,:)=cell2mat(tbl_delay(2,6));
    num_trial(i)=length(all_cuedelay)/2;
    rank_cuerate(i,:)=sort(temp_cuerate/max(temp_cuerate),'descend');
    rank_samplerate(i,:)=sort(temp_samplerate/max(temp_samplerate),'descend');
    rank_delayrate(i,:)=sort(temp_delayrate/max(temp_delayrate),'descend');
  %  store_cuerate(i)=mean(temp_cuerate); %mean firing rate for each cell
  %  store_samplerate(i)=mean(temp_samplerate);
    store_cuesi(i)=max(temp_cuerate)/((sum(temp_cuerate)-max(temp_cuerate))/(length(location_to_fit)-1));  %selectivity index for each cell
    store_samplesi(i)=max(temp_samplerate)/((sum(temp_samplerate)-max(temp_samplerate))/(length(location_to_fit)-1));
    store_delaysi(i)=max(temp_delayrate)/((sum(temp_delayrate)-max(temp_delayrate))/(length(location_to_fit)-1));
  %  if num_trial(i)>90
  %     good_data(i)=0;  
  %  end
    else
    good_data(i)=0; 
  %  stim_p(i)=[];
  %  delay_p(i)=[];
    end

end
    cue_p(cue_p==0) = NaN;
    sample_p(sample_p==0) = NaN;
    cuevsdelay_p(cuevsdelay_p==0) = NaN;
    delay_p(delay_p==0)=NaN;
    
    data_inclusion_index=find(good_data==1);
    find_bothmain=intersect(intersect(find(cuevsdelay_p(:,1)<0.05),find(cuevsdelay_p(:,2)<0.05)),data_inclusion_index);   
    find_nomain=intersect(intersect(find(cuevsdelay_p(:,1)>0.05),find(cuevsdelay_p(:,2)>0.05)),data_inclusion_index);
    find_sigstim_main1=intersect(find(cuevsdelay_p(:,1)<0.05),data_inclusion_index);
    find_sigstim_main2=intersect(find(cuevsdelay_p(:,2)<0.05),data_inclusion_index);
    find_sigstim_interaction=intersect(find(cuevsdelay_p(:,3)<0.05),data_inclusion_index);
   
    
    main1_resample_num=floor(resample_size*Parietal_proportion_resample(1));
    main2_resample_num=floor(resample_size*Parietal_proportion_resample(2));
    for r=1:500  %permutation test loop
    resample_main1=find_sigstim_main1(randperm(length(find_sigstim_main1),main1_resample_num));
    main2_included=intersect(find_sigstim_main2,resample_main1);
    main2_left=setdiff(find_sigstim_main2,find_sigstim_main1);
    if length(main2_included)>=main2_resample_num
        resample_main2=main2_included(randperm(length(main2_included),main2_resample_num));
    else
        resample_main2=[main2_included;main2_left(randperm(length(main2_left),main2_resample_num-length(main2_included)))];
    end
    resample_all_part1=unique([resample_main1;resample_main2]);
    resample_all_part2=find_nomain(randperm(length(find_nomain),resample_size-length(resample_all_part1)));
    resample_all=[resample_all_part1;resample_all_part2];
    interaction_inresample_parietalproportion(r)=length(intersect(resample_all,find_sigstim_interaction))/resample_size;
    end
    
    main1_resample_num=floor(resample_size*PFC_proportion_resample(1));
    main2_resample_num=floor(resample_size*PFC_proportion_resample(2));
    for r=1:500
    resample_main1=find_sigstim_main1(randperm(length(find_sigstim_main1),main1_resample_num));
    main2_included=intersect(find_sigstim_main2,resample_main1);
    main2_left=setdiff(find_sigstim_main2,find_sigstim_main1);
    if length(main2_included)>=main2_resample_num
        resample_main2=main2_included(randperm(length(main2_included),main2_resample_num));
    else
        resample_main2=[main2_included;main2_left(randperm(length(main2_left),main2_resample_num-length(main2_included)))];
    end
    resample_all_part1=unique([resample_main1;resample_main2]);
    resample_all_part2=find_nomain(randperm(length(find_nomain),resample_size-length(resample_all_part1)));
    resample_all=[resample_all_part1;resample_all_part2];
    interaction_inresample_pfcproportion(r)=length(intersect(resample_all,find_sigstim_interaction))/resample_size;
    end

    figure;
    histogram(interaction_inresample_pfcproportion,'BinWidth',0.01);
    hold on;
    histogram(interaction_inresample_parietalproportion,'BinWidth',0.01);
    legend('pfc proportion','parietal proportion');
    
    
    find_sigstim_CS1=setdiff(setdiff(find_sigstim_main1,find_sigstim_main2),find_sigstim_interaction);
    find_sigstim_CS2=setdiff(setdiff(find_sigstim_main2,find_sigstim_main1),find_sigstim_interaction);
    find_sigstim_informative=unique([find_sigstim_main1;find_sigstim_main2;find_sigstim_interaction]);
    find_sigstim_LMS=setdiff(intersect(find_sigstim_main1,find_sigstim_main2),find_sigstim_interaction);
    find_sigstim_NMS=find_sigstim_interaction;
    find_sigstim_NMS1=intersect(find_sigstim_interaction,intersect(find_sigstim_main1,find_sigstim_main2));
    find_sigstim_NMS2=intersect(find_sigstim_interaction,setdiff(find_sigstim_main1,find_sigstim_main2));
    find_sigstim_NMS3=intersect(find_sigstim_interaction,setdiff(find_sigstim_main2,find_sigstim_main1));
    find_sigstim_NMS4=setdiff(find_sigstim_interaction,unique([find_sigstim_main1;find_sigstim_main2]));
    find_sigstim_CS=[find_sigstim_CS1;find_sigstim_CS2];
      

    prop_sigstim_main1=length(find_sigstim_main1)/length(data_inclusion_index);
    prop_sigstim_main2=length(find_sigstim_main2)/length(data_inclusion_index);
    prop_sigstim_interaction=length(find_sigstim_interaction)/length(data_inclusion_index);
    
    main_effect_y=[prop_sigstim_main1,prop_sigstim_main2,prop_sigstim_interaction];
    f=figure;
    bar(main_effect_y);
    %h.EdgeColor = 'none';
    %h.FaceColor=[];
    hold on;
    xlimit = get(gca,'xlim');
    plot(xlimit,[0.05 0.05]);
  %  title('Prefrontal location vs epoch');
    set(gca,'xticklabel',{'location','cue vs delay','interaction'})
   % legend('stim','delay');
    ylim([0,0.5]);
  %  text(1,0.4,strcat('n=',num2str(length(data_inclusion_index))));
  yticks([0:0.1:0.5])
  yticklabels({'0','10','20','30','40','50'});  
  set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
  f.Position(3:4)=[600,400]; %figure size
  box off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%line plot for p value dring cue and
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%delay period
    count_category=[0,0,0,0];
    figure;   %plot f-score/p-value change
    hold on
    category1_index=[];
    category2_index=[];
    category3_index=[];
    category4_index=[];
    category1_teststat=[];
    category2_teststat=[];
    category3_teststat=[];
    category4_teststat=[];
    for i=1:length(data_inclusion_index)
        if cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))>=0.05
            plot1=plot([f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i))],'-ok');
          %  plot1=plot([cue_p(data_inclusion_index(i)),sample_p(data_inclusion_index(i))],'-ok');
            plot1.Color(4) = 0.05;
            count_category(1)=count_category(1)+1;
            category1_teststat=[category1_teststat;[cue_p(data_inclusion_index(i)),delay_p(data_inclusion_index(i))]];
            category1_index=[category1_index,data_inclusion_index(i)];
        elseif cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))>=0.05
            plot2=plot([f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i))],'--sr');
          %  plot2=plot([cue_p(data_inclusion_index(i)),sample_p(data_inclusion_index(i))],'--sr');
            count_category(2)=count_category(2)+1;
            category2_teststat=[category2_teststat;[cue_p(data_inclusion_index(i)),delay_p(data_inclusion_index(i))]];
            category2_index=[category2_index,data_inclusion_index(i)];
        elseif cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))<0.05
            plot3=plot([f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i))],'--dm');
          %  plot3=plot([cue_p(data_inclusion_index(i)),sample_p(data_inclusion_index(i))],'--dm');
            count_category(3)=count_category(3)+1;
            category3_teststat=[category3_teststat;[cue_p(data_inclusion_index(i)),delay_p(data_inclusion_index(i))]];
            category3_index=[category3_index,data_inclusion_index(i)];
        else cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))<0.05
            plot4=plot([f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i))],':xg'); 
          %  plot4=plot([cue_p(data_inclusion_index(i)),sample_p(data_inclusion_index(i))],':xg');
            count_category(4)=count_category(4)+1;
            category4_teststat=[category4_teststat;[cue_p(data_inclusion_index(i)),delay_p(data_inclusion_index(i))]];
            category4_index=[category4_index,data_inclusion_index(i)];
        end    
    end
  %%%%%%%%%%%%%%%%select cells to plot psth to show example
    [B,best_index]=sort(category4_teststat(:,2),'ascend');
    cell_index=select_index(category4_index(best_index(1:10)));
  %  category_save=select_index(index_category4); %save certain category for certain area
  %  save prefrontal_category4.mat category_save
  category_save{1}=select_index(find_sigstim_CS1); 
  category_save{2}=select_index(find_sigstim_CS2); 
  category_save{3}=select_index(find_sigstim_LMS);
  category_save{4}=select_index(find_sigstim_NMS); 
 % save parietal_category.mat category_save
  proportion_category=count_category/length(data_inclusion_index);
    legend([plot1,plot2,plot3,plot4],{['non-selective',' ',num2str(proportion_category(1)*100),'%.2f'],...
        ['selectivity disappear',' ',num2str(proportion_category(2)*100),'%.2f'],...
        ['selectivity emerge',' ',num2str(proportion_category(3)*100),'%.2f'],...
        ['duo-selective',' ',num2str(proportion_category(4)*100),'%.2f']});  
    title('Parietal f value for location');
    set(gca,'XTick',[]);
    xticks([1 2]);
 %   xticklabels({'CUE','SAMPLE'});
    ylim([0,30]);
 
 
 
 
  %scatter plots for f value for location selectivity during cue and delay
  %only color plot the NMS cells/all cells
  temp_count=[0,0,0,0];
  f=figure;   %plot f-score/p-value change
  hold on
    for i=1:length(data_inclusion_index)
      %  if ismember(data_inclusion_index(i),find_sigstim_NMS) & cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))<0.05
        if cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))<0.05      
            plot4=scatter(f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i)),100,'xg','LineWidth',2); 
            temp_count(4)=temp_count(4)+1;
        %elseif ismember(data_inclusion_index(i),find_sigstim_NMS) & cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))>=0.05
        elseif cue_p(data_inclusion_index(i))<0.05 & delay_p(data_inclusion_index(i))>=0.05       
            plot2=scatter(f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i)),100,'xr','LineWidth',2);
            temp_count(2)=temp_count(2)+1;
       % elseif ismember(data_inclusion_index(i),find_sigstim_NMS) & cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))<0.05
        elseif cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))<0.05
            plot3=scatter(f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i)),100,'xm','LineWidth',2);
            temp_count(3)=temp_count(3)+1;
       % elseif ismember(data_inclusion_index(i),find_sigstim_NMS) & cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))>=0.05
        elseif cue_p(data_inclusion_index(i))>=0.05 & delay_p(data_inclusion_index(i))>=0.05
            plot1=scatter(f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i)),20,'ok','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'LineWidth',2);
            temp_count(1)=temp_count(1)+1;
        else 
            plot5=scatter(f_cue(data_inclusion_index(i)),f_delay(data_inclusion_index(i)),20,'ok','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'LineWidth',2);
        end    
    end
NMS_proportion_category(1)=temp_count(1)/sum(temp_count);
NMS_proportion_category(2)=temp_count(2)/sum(temp_count);
NMS_proportion_category(3)=temp_count(3)/sum(temp_count);
NMS_proportion_category(4)=temp_count(4)/sum(temp_count);
lgd=legend([plot1,plot2,plot3,plot4],{['non-selective',' ',num2str(round(NMS_proportion_category(1)*1000)/10),'%'],...
        ['cue selective',' ',num2str(round(NMS_proportion_category(2)*1000)/10),'%'],...
        ['delay selective',' ',num2str(round(NMS_proportion_category(3)*1000)/10),'%'],...
        ['duo-selective',' ',num2str(round(NMS_proportion_category(4)*1000)/10),'%']}); 
lgd.FontSize=18;
legend boxoff;   
set(gca,'fontsize',16,'FontWeight','bold','LineWidth',2); %label size
f.Position(3:4)=[560,420]; %figure size
%box off;
ylim([0,15]);
xlim([0,15]);


 %%%%%%%%%%%%%%%%%%%same as figure3, but just use NMS cells   
    figure; 
    hold on;
    for i=1:length(find_sigstim_NMS) %overlay NMS cells
        if cue_p(find_sigstim_NMS(i))>=0.05 & delay_p(find_sigstim_NMS(i))>=0.05
            plot5=plot([f_cue(find_sigstim_NMS(i)),f_delay(find_sigstim_NMS(i))],'-ok');
         %  plot1=plot([cue_p(find_sigstim_NMS(i)),sample_p(find_sigstim_NMS(i))],'-ok');  
            plot5.Color(4) = 0.15;
        elseif cue_p(find_sigstim_NMS(i))<0.05 & delay_p(find_sigstim_NMS(i))>=0.05
            plot6=plot([f_cue(find_sigstim_NMS(i)),f_delay(find_sigstim_NMS(i))],'--sr');
         %  plot2=plot([cue_p(find_sigstim_NMS(i)),sample_p(find_sigstim_NMS(i))],'--sr');
        elseif cue_p(find_sigstim_NMS(i))>=0.05 & delay_p(find_sigstim_NMS(i))<0.05
            plot7=plot([f_cue(find_sigstim_NMS(i)),f_delay(find_sigstim_NMS(i))],'--dm');
         %  plot3=plot([cue_p(find_sigstim_NMS(i)),sample_p(find_sigstim_NMS(i))],'--dm');
        else cue_p(find_sigstim_NMS(i))<0.05 & delay_p(find_sigstim_NMS(i))<0.05
            plot8=plot([f_cue(find_sigstim_NMS(i)),f_delay(find_sigstim_NMS(i))],':xg');
         %  plot4=plot([cue_p(find_sigstim_NMS(i)),sample_p(find_sigstim_NMS(i))],':xg');
        end
    end
   % ylim([0,30]);
    NMS_proportion_category(1)=length(intersect(find_sigstim_NMS,index_category1))/length(find_sigstim_NMS);
    NMS_proportion_category(2)=length(intersect(find_sigstim_NMS,index_category2))/length(find_sigstim_NMS);
    NMS_proportion_category(3)=length(intersect(find_sigstim_NMS,index_category3))/length(find_sigstim_NMS);
    NMS_proportion_category(4)=length(intersect(find_sigstim_NMS,index_category4))/length(find_sigstim_NMS);
    legend([plot5,plot6,plot7,plot8],{['non-selective',' ',num2str(NMS_proportion_category(1)*100),'%'],...
        ['selectivity disappear',' ',num2str(NMS_proportion_category(2)*100),'%'],...
        ['selectivity emerge',' ',num2str(NMS_proportion_category(3)*100),'%'],...
        ['duo-selective',' ',num2str(NMS_proportion_category(4)*100),'%']}); 
    title('Parietal f value for location, NMS only');
    set(gca,'XTick',[]);
    xticks([1 2]);
    xticklabels({'CUE','SAMPLE'});
    
    %%%%%%%%%%%plot prefered location change for cells that are selectivie
    %%%%%%%%%%%in both cue and delay period
    figure; 
    hold on
    c=1;
    for i=1:length(good_data)

        if cue_p(i)<0.05 & delay_p(i)<0.05 & gof_cue_r2(i)>0.3 & gof_delay_r2(i)>0.3
       %     scatter(gauss_cue(i,2),gauss_delay(i,2),'ok');
       % elseif cue_p(i)<0.05 & delay_p(i)>=0.05
       %     scatter(gauss_cue(i,2),gauss_sample(i,2),'sr');
       % elseif cue_p(i)>=0.05 & delay_p(i)<0.05
       %     scatter(gauss_cue(i,2),gauss_delay(i,2),'dm');
       % else
           scatter(gauss_cue(i,2),gauss_delay(i,2),60,[0 0.447 0.741],'LineWidth',1.5); 
          % lsline;
           store_gauss_cue(c)=gauss_cue(i,2);
           store_gauss_delay(c)=gauss_delay(i,2);
           c=c+1;
            %scatter(raw_cue_preferloc(i),raw_sample_preferloc(i),80,'xg'); 
        end    
    end
    xlim([-2,2]);
    ylim([-2,2]);
   % title('Parietal duo-selective prefered location');
    %xlabel('cue prefer location');
    xticklabels({'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'});
    %ylabel('sample prefer location');
    yticklabels({'-90','-67.5','-45','-22.5','0','22.5','45','67.5','90'});
    hold on;
    plot([-2,2],[-2,2],'k','LineWidth',1);
    g=sum((store_gauss_cue-store_gauss_delay).^2);
    y=sum((store_gauss_delay-mean(store_gauss_delay)).^2); 
    r2=1-g/y;
    hAx=gca; 
    hAx.LineWidth=1.5;hAx.FontSize = 12;
    t = text(0.7,-1.5,'R squared=0.761');
    t.FontSize = 12;
    save PFC_gauss_cue.mat store_gauss_cue
    save PFC_gauss_delay.mat store_gauss_delay
    
    figure; % plot prefered location for new generated selective cells
    new_selective=find(cue_p>=0.05 & delay_p<0.05);
    %histogram(gauss_sample(new_selective,2),[-2:0.2:2],'FaceColor',[1,0,1]);
    histogram(raw_delay_preferloc(new_selective),[0:10],'FaceColor',[1,0,1]);
    title('Parietal emerge-selective prefered location');
    ylabel('count');
    xlabel('preferred location in sample period');
    xticks([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]);
    xticklabels({'-90','-45','-22.5','-11.25','0','11.25','22.5','45','90'});

    %plot prefered location for cue and delay period respectively
    figure; 
    new_selective=find(cue_p<0.05);
    histogram(raw_cue_preferloc(new_selective),[0:10],'FaceAlpha',0.2);
    hold on;
    new_selective=find(delay_p<0.05);
    histogram(raw_delay_preferloc(new_selective),[0:10],'FaceAlpha',0.2);
    legend('cue preferred location','sample preferred location');
    xlabel('preferred location');
    xticks([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]);
    xticklabels({'-90','-45','-22.5','-11.25','0','11.25','22.5','45','90'});
    title('Parietal all prefered location');
    %{
    figure; %quantify tuning sharpness by normalize ranking
    new_selective=find(cue_p<0.05 & sample_p<0.05);
    plot(mean(rank_cuerate(new_selective,:)));
    hold on;
    plot(mean(rank_samplerate(new_selective,:)));
    title('Parietal rank-ordered response');
    legend('CUE','SAMPLE');
    %}
