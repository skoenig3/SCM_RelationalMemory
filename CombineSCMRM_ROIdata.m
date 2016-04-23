function CombineSCMRM_ROIdata(cortexfiles,data_dir,figure_dir)
%written by Seth Konig 9/8/15
%Code combines ROI data for the SCMRM task.


if strcmp(cortexfiles{1},'All') %if analyzing data across monkeys.
    num_points = NaN(length(cortexfiles)-1,4);
    all_sets_numFix = NaN(length(cortexfiles)-1,4);
    all_sets_10numFix = NaN(length(cortexfiles)-1,4);
    all_sets_Time = NaN(length(cortexfiles)-1,4);
    all_sets_Time3 = NaN(length(cortexfiles)-1,4);
    all_sets_Time15 = NaN(length(cortexfiles)-1,4);
    all_time_sets = zeros(5,7000);
    
    numsets = length(cortexfiles)-1;
    
    for sets = 2:length(cortexfiles)
        load([data_dir cortexfiles{sets}(1:8) '_' cortexfiles{sets}(10) '-ROIdata.mat'])
        
        for type = 1:4
            %nans appear if monkey looked way during the period of interest
            %(e.g. 500-1500 ms)
            all_sets_numFix(sets-1,type) = 100*nanmean(allnumFixationsInROI{type});
            all_sets_10numFix(sets-1,type) = 100*nanmean(num10FixationsInROI{type});
            all_sets_Time(sets-1,type) = 100*nanmean(allTime_ROI_timeWindow{type});
            all_sets_Time3(sets-1,type) = 100*nanmean(Time3_ROI_timeWindow{type});
            all_sets_Time15(sets-1,type) = 100*nanmean(Time15_ROI_timeWindow{type});
        end
        all_time_sets = all_time_sets +all_time;
    end
    all_time_sets = all_time_sets/numsets;
    
    
    p_numFix = NaN(1,4);
    p_numFix10 = NaN(1,4);
    p_Time = NaN(1,4);
    p_Time3 = NaN(1,4);
    p_Time15 = NaN(1,4);
    for type = [1 3 4]
        [~,p_numFix(type)] = ttest2(all_sets_numFix(:,type),all_sets_numFix(:,2));
        [~,p_numFix10(type)] = ttest2(all_sets_10numFix(:,type),all_sets_10numFix(:,2));
        [~,p_Time(type)] = ttest2(all_sets_Time(:,type),all_sets_Time(:,2));
        [~,p_Time3(type)] = ttest2(all_sets_Time3(:,type),all_sets_Time3(:,2));
        [~,p_Time15(type)] = ttest2(all_sets_Time15(:,type),all_sets_Time15(:,2));
    end
    
    
    figure
    subplot(2,3,1)
    hold on
    bar(mean(all_sets_10numFix));
    errorb(mean(all_sets_10numFix),std(all_sets_10numFix)./sqrt(numsets))
    for type = [1 3 4]
        if p_numFix10(type) < 0.05
            plot(type,mean(all_sets_10numFix(:,type))+std(all_sets_10numFix(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of 1st 10 Fixations in ROI')
    
    subplot(2,3,2)
    hold on
    bar(mean(all_sets_Time3));
    errorb(mean(all_sets_Time3),std(all_sets_Time3)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time3(type) < 0.05
            plot(type,mean(all_sets_Time3(:,type))+std(all_sets_Time3(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of Time in 1st 3 seconds in ROI')
    
    subplot(2,3,3)
    hold on
    bar(mean(all_sets_Time15));
    errorb(mean(all_sets_Time15),std(all_sets_Time15)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time15(type) < 0.05
            plot(type,mean(all_sets_Time15(:,type))+std(all_sets_Time15(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of Time in 0.5-1.5 sec window in ROI')
    
    subplot(2,3,4)
    hold on
    bar(mean(all_sets_numFix));
    errorb(mean(all_sets_numFix),std(all_sets_numFix)./sqrt(numsets))
    for type = [1 3 4]
        if p_numFix(type) < 0.05
            plot(type,mean(all_sets_numFix(:,type))+std(all_sets_numFix(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of All Fixations in ROI')
    
    subplot(2,3,5)
    hold on
    bar(mean(all_sets_Time));
    errorb(mean(all_sets_Time),std(all_sets_Time)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time(type) < 0.05
            plot(type,mean(all_sets_Time(:,type))+std(all_sets_Time(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of All Time in ROI')
    
    clrs = 'brmgg';
    line = {'-','-','-','-','--'};
    subplot(2,3,6)
    hold on
    for r = 1:size(all_time,1)
        plot(100*filtfilt(1/200*ones(1,200),1,all_time_sets(r,:)),...
            [clrs(r) line{r}]);
    end
    hold off
    xlabel('Time (ms)')
    ylabel('% of time in ROI')
    legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
    title('% of Time in ROI over Time')
    yl = ylim;
    ylim([0 yl(2)])
    
    save_and_close_fig(figure_dir,['Allmonkeys_allsets'])
    
else
    num_points = NaN(length(cortexfiles),4);
    all_sets_numFix = NaN(length(cortexfiles),4);
    all_sets_10numFix = NaN(length(cortexfiles),4);
    all_sets_Time = NaN(length(cortexfiles),4);
    all_sets_Time3 = NaN(length(cortexfiles),4);
    all_sets_Time15 = NaN(length(cortexfiles),4);
    all_time_sets = zeros(5,7000);
    
    numsets = length(cortexfiles);
    
    for sets = 1:length(cortexfiles)
        load([data_dir cortexfiles{sets}(1:8) '_' cortexfiles{sets}(10) '-ROIdata.mat'])
        
        for type = 1:4
            num_points(sets,type)= sum(~isnan(allnumFixationsInROI{type}));
            all_sets_numFix(sets,type) = 100*nanmean(allnumFixationsInROI{type});
            all_sets_10numFix(sets,type) = 100*nanmean(num10FixationsInROI{type});
            all_sets_Time(sets,type) = 100*nanmean(allTime_ROI_timeWindow{type});
            all_sets_Time3(sets,type) = 100*nanmean(Time3_ROI_timeWindow{type});
            all_sets_Time15(sets,type) = 100*nanmean(Time15_ROI_timeWindow{type});
        end
        all_time_sets = all_time_sets +all_time;
    end
    all_time_sets = all_time_sets/numsets;
    
    
    p_numFix = NaN(1,4);
    p_numFix10 = NaN(1,4);
    p_Time = NaN(1,4);
    p_Time3 = NaN(1,4);
    p_Time15 = NaN(1,4);
    for type = [1 3 4]
        [~,p_numFix(type)] = ttest2(all_sets_numFix(:,type),all_sets_numFix(:,2));
        [~,p_numFix10(type)] = ttest2(all_sets_10numFix(:,type),all_sets_10numFix(:,2));
        [~,p_Time(type)] = ttest2(all_sets_Time(:,type),all_sets_Time(:,2));
        [~,p_Time3(type)] = ttest2(all_sets_Time3(:,type),all_sets_Time3(:,2));
        [~,p_Time15(type)] = ttest2(all_sets_Time15(:,type),all_sets_Time15(:,2));
    end
    
    
    figure
    subplot(2,3,1)
    hold on
    bar(mean(all_sets_10numFix));
    errorb(mean(all_sets_10numFix),std(all_sets_10numFix)./sqrt(numsets))
    for type = [1 3 4]
        if p_numFix10(type) < 0.05
            plot(type,mean(all_sets_10numFix(:,type))+std(all_sets_10numFix(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of 1st 10 Fixations in ROI')
    
    subplot(2,3,2)
    hold on
    bar(mean(all_sets_Time3));
    errorb(mean(all_sets_Time3),std(all_sets_Time3)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time3(type) < 0.05
            plot(type,mean(all_sets_Time3(:,type))+std(all_sets_Time3(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of Time in 1st 3 seconds in ROI')
    
    subplot(2,3,3)
    hold on
    bar(mean(all_sets_Time15));
    errorb(mean(all_sets_Time15),std(all_sets_Time15)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time15(type) < 0.05
            plot(type,mean(all_sets_Time15(:,type))+std(all_sets_Time15(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of Time in 0.5-1.5 sec window in ROI')
    
    subplot(2,3,4)
    hold on
    bar(mean(all_sets_numFix));
    errorb(mean(all_sets_numFix),std(all_sets_numFix)./sqrt(numsets))
    for type = [1 3 4]
        if p_numFix(type) < 0.05
            plot(type,mean(all_sets_numFix(:,type))+std(all_sets_numFix(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of All Fixations in ROI')
    
    subplot(2,3,5)
    hold on
    bar(mean(all_sets_Time));
    errorb(mean(all_sets_Time),std(all_sets_Time)./sqrt(numsets))
    for type = [1 3 4]
        if p_Time(type) < 0.05
            plot(type,mean(all_sets_Time(:,type))+std(all_sets_Time(:,type)),'k*')
        end
    end
    hold off
    set(gca,'Xtick',1:4)
    set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
    ylabel('Percentage')
    title('% of All Time in ROI')
    
    clrs = 'brmgg';
    line = {'-','-','-','-','--'};
    subplot(2,3,6)
    hold on
    for r = 1:size(all_time,1)
        plot(100*filtfilt(1/200*ones(1,200),1,all_time_sets(r,:)),...
            [clrs(r) line{r}]);
    end
    hold off
    xlabel('Time (ms)')
    ylabel('% of time in ROI')
    legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
    title('% of Time in ROI over Time')
    yl = ylim;
    ylim([0 yl(2)])
    
    save_and_close_fig(figure_dir,[cortexfiles{1}(1:2) '_allsets'])
    
end

