for SS =12:length(preNIN_pha1_2)
    LFP_SS       = importdata(['RAWDATA/RawLFP/Phase1/',preNIN_pha1_2{SS}]);
    for TRIAL=1:length(LFP_SS)
        temp = LFP_SS{TRIAL};
        if isfield(temp.channels(1),'lfp_filter')
            if (length(temp.channels)==64)
                for ch=1:(length(temp.channels) -1)
                    TEM(ch,:) = temp.channels(ch).lfp_filter;
                end
            end
        end
        TEM1               = zeros(size(TEM,1),size(TEM,1));
        for i=1:size(TEM,1)
            y             = TEM(i,:)';
            if (i==1)
                X             = TEM(2:size(TEM,1),:)';
                [b,dev,stats] = glmfit(X,y);
                TEM1(i,2:end)  = b(2:end,:)';
                clear b
            elseif (i==size(TEM,1))
                X             = TEM(1:(size(TEM,1)-1),:)';
                [b,dev,stats] = glmfit(X,y);
                TEM1(i,1:size(TEM,1)-1)   = b(2:size(TEM,1),:)';
                clear b
            else
                a             = [1:i-1 i+1:size(TEM,1)];
                X             = TEM(a,:)';
                [b,dev,stats] = glmfit(X,y);
                TEM1(i,a)      = b(2:size(TEM,1),:)';
                clear b
            end
        end
        % % %% Indentifying the noise channel based on collinearity
        for i=1:size(TEM1,1)
            if(-1<TEM1(i,:))
                TTC(i,:)      = 1;
            elseif(TEM1(i,:)>1)
                TTC(i,:)      = 1;
            else
                TTC(i,:)      = 0;
            end
        end
        % % %% Replace NAN for noisy channel
        for ch=1:size(TEM1,1)
            if(TTC(ch,:)==0)
                LFP(ch,:) = zeros(1,size(TEM,2));
            else
                LFP(ch,:)= TEM(ch,:);
            end
        end
        save(['FilterData/Phase1/Active/NIN/',preNIN_pha1_2{SS}(1:8),'_',num2str(TRIAL),'.mat'],'LFP','-v7.3')
        clear temp TEM1 TTC LFP TEM;
    end
    clear LFP_SS;
end