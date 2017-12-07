clar
files = {'TO150901.2','TO150903.2','TO150904.3',...
    'TO150908.2','TO150909.2','TO150910.2','TO150911.2',...
    'TO150915.2','TO150916.2','TO150918.2',...
    'TO151001.2','TO151002.2','TO151006.2','TO151007.2',...
    'PW150901.2','PW150902.2','PW150903.2','PW150904.2',...
    'PW150908.2','PW150910.2','PW150911.2',...
    'PW150914.2','PW150916.2','PW150917.2','PW150918.2',...
    'PW150921.2','PW150922.2','PW150923.2','PW150924.2'};

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';%location of processed eye data
ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\'; %locations of regions of interest
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Figures\'; %where to put plots
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\'; %where all the image sets are

imageX = 800;
imageY = 600;

dist_threshold = 48;
N_fix = 40;

all_img_types = NaN(1,length(files)*36);
all_recurrence_rates = NaN(1,length(files)*36);
all_laminarity = NaN(1,length(files)*36);
all_forward_trace = NaN(1,length(files)*36);
all_reverse_trace = NaN(1,length(files)*36);
all_ROI_time = NaN(1,length(files)*36);
all_ROI_515 = NaN(1,length(files)*36);
all_ROI_first = NaN(1,length(files)*36);
mean_fix_durs = NaN(1,length(files)*36);
all_corm = NaN(1,length(files)*36);

img_id = 0;
for sets = 1:length(files)
    
    
    %---Eye Movement Analysis---%
    load([data_dir files{sets}(1:8) '_' files{sets}(10) '-fixation.mat'])
    
    %---ROI anlalysis---%
    load([ROI_dir itemfile(1:end-4) '_ROIs.mat']);
    
    setnum = str2double(itemfile(6:7));
    
    trial_types = NaN(1,72);
    fixdurs = NaN(72,50); %each trials fixation durations
    sacamps = NaN(72,50); %each trials saccade amplitudes
    
    for img = 1:36;
        if img == 18 && setnum == 30%image accidentaly was shown in a previous set
            continue
        end
        
        %---Grab Important Vairables---%
        %for novel images
        nov_allval = per(2*img-1).allval;
        nov_alltim = per(2*img-1).alltim;
        nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
        
        if nov_img_off-nov_img_on > 14000
            continue
        end
        
        nov_x = fixationstats{2*img-1}.XY(1,nov_img_on:nov_img_off);
        nov_y = fixationstats{2*img-1}.XY(2,nov_img_on:nov_img_off);
        nov_pupil = pupildata{2*img-1}(round(nov_img_on/5):round(nov_img_off/5));
        
        nov_fix = fixationstats{2*img-1}.fixations;
        nov_fixtimes = fixationstats{2*img-1}.fixationtimes;
        nov_sactimes = fixationstats{2*img-1}.saccadetimes;
        
        pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
        post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
        pre_img_saccades = find(nov_sactimes(1,:) <= nov_img_on);
        post_img_saccades = find(nov_sactimes(2,:) > nov_img_off);
        
        nov_fix(:,post_img_fix) = [];
        nov_fix(:,pre_img_fix) = [];
        nov_fixtimes(:,post_img_fix) = [];
        nov_fixtimes(:,pre_img_fix ) = [];
        nov_sactimes(:,post_img_saccades) = [];
        nov_sactimes(:,pre_img_saccades) = [];
        
        
        %for repeat images
        rep_allval = per(2*img).allval;
        rep_alltim = per(2*img).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        if rep_img_off-rep_img_on > 14000
            continue
        end
        
        rep_x = fixationstats{2*img}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{2*img}.XY(2,rep_img_on:rep_img_off);
        
        rep_fix = fixationstats{2*img}.fixations;
        rep_fixtimes = fixationstats{2*img}.fixationtimes;
        rep_sactimes = fixationstats{2*img}.saccadetimes;
        
        pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
        post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
        pre_img_saccades = find(rep_sactimes(1,:) <= rep_img_on);
        post_img_saccades = find(rep_sactimes(2,:) > rep_img_off);
        
        rep_fix(:,post_img_fix) = [];
        rep_fix(:,pre_img_fix) = [];
        rep_fixtimes(:,post_img_fix) = [];
        rep_fixtimes(:,pre_img_fix) = [];
        rep_sactimes(:,post_img_saccades) = [];
        rep_sactimes(:,pre_img_saccades) = [];
        
        %%
        rep_trialtype = trialtype(2,img);
        %---Import ROIs---%
        if rep_trialtype == 2 || rep_trialtype == 3 %familiar and replaced, respectively
            ROI1 = ROIs{img}{1};
            ROI1(1) = ROI1(1)-36;
            ROI1(2) = ROI1(2)+36;
            ROI1(3) = ROI1(3)-36;
            ROI1(4) = ROI1(4)+36;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 =ROI1;
        elseif rep_trialtype == 4 %moved
            ROI1 = ROIs{img}{1};%original location
            ROI1(1) = ROI1(1)-36;
            ROI1(2) = ROI1(2)+36;
            ROI1(3) = ROI1(3)-36;
            ROI1(4) = ROI1(4)+36;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 = ROIs{img}{2};%new location after move
            ROI2(1) = ROI2(1)-36;
            ROI2(2) = ROI2(2)+36;
            ROI2(3) = ROI2(3)-36;
            ROI2(4) = ROI2(4)+36;
            ROI2(ROI2 < 1) = 1;
            ROI2(ROI2 > imageX) = imageX;
            if ROI2(3) > imageY; ROI2(3) = imageY;end
            if ROI2(4) > imageY; ROI2(4) = imageY;end
        end
        
        %---Find Fixations within ROI---%
        ROIfix1 = find(...
            (nov_fix(1,:) > ROI1(1) & nov_fix(1,:) < ROI1(2)) & ...
            (nov_fix(2,:) > ROI1(3) & nov_fix(2,:) < ROI1(4))); %ROIfix1 contains indices of nov_fix corresponding to fixations that occur in ROI of novel image.
        ROIfix2 = find(...
            rep_fix(1,:) > ROI2(1) & rep_fix(1,:) < ROI2(2) & ...
            rep_fix(2,:) > ROI2(3) & rep_fix(2,:) < ROI2(4));%ROIfix2 contains indices of rep_fix corresponding to fixations that occur in ROI of the manipulated image
        
        if rep_trialtype ==4 %moved image. Data from the old ROI (where the object was before being moved)
            ROIfix2_oldROI = find(...
                (rep_fix(1,:) > ROI1(1) & rep_fix(1,:) < ROI1(2) & ...
                rep_fix(2,:) > ROI1(3) & rep_fix(2,:) < ROI1(4)));%ROIfix2_oldROI contains indices of rep_fix corresponding to fixations that occur in the 'old' ROI of the 'moved' image
        else
            ROIfix2_oldROI = [];
        end
        
        if ~isempty(ROIfix1)
            %---Calculate Area of ROI---%
            %only save area of data points so structure holds
            area_ROI1 = abs(ROI1(1)-ROI1(2))*abs(ROI1(3)-ROI1(4)); %area of ROI for novel image
            area_ROI2 = abs(ROI2(1)-ROI2(2))*abs(ROI2(3)-ROI2(4)); %for replaced, moved or familiar
            
            %---Calculate Average Time in ROI---%
            tempvec1=zeros(1,7000);%novel image
            for j=1:length(ROIfix1)
                tempvec1(nov_fixtimes(1,ROIfix1(j)):nov_fixtimes(2,ROIfix1(j)))=1;
            end
            
            tempvec2=zeros(1,7000);%manipulated image new ROI
            for j=1:length(ROIfix2)
                tempvec2(rep_fixtimes(1,ROIfix2(j)):rep_fixtimes(2,ROIfix2(j)))=1;
            end
            
            tempvec3=zeros(1,7000);%manipulated imaged old moved ROI
            for j=1:length(ROIfix2_oldROI)
                tempvec3(rep_fixtimes(1,ROIfix2_oldROI(j)):rep_fixtimes(2,ROIfix2_oldROI(j)))=1;
            end
            
            %---Calculate Stats for the Novel image---%
            nov_propall_fixations = length(ROIfix1)/length(nov_fix);
            nov_prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
            nov_propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
            nov_prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
            nov_prop15_time = nansum(tempvec1(500:1500))/sum(~isnan(tempvec1(500:1500)));
            
            %---Calculate Stats for the Second Presentation---%
            if trialtype == 4 %for moved
                rep_propall_fixations = (length(ROIfix2)+length(ROIfix2_oldROI))/length(rep_fix);
                rep_prop10_fixaitons = (length(ROIfix2(ROIfix2 <= 10))+...
                    length(ROIfix2_oldROI(ROIfix2_oldROI <= 10)))./10;
                rep_propall_time = (nansum(tempvec2)+nansum(tempvec3))/sum(~isnan(tempvec2));
                rep_prop3_time = (nansum(tempvec2(1:3000))+nansum(tempvec3(1:3000)))...
                    /sum(~isnan(tempvec2(1:3000)));
                rep_prop15_time = (nansum(tempvec2(500:1500))+nansum(tempvec3(500:1500)))...
                    /sum(~isnan(tempvec2(500:1500)));
                
            else %for replaced or repeat
                rep_propall_fixations = length(ROIfix2)/length(rep_fix);
                rep_prop10_fixaitons = length(ROIfix2(ROIfix2 <= 10))/10;
                rep_propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
                rep_prop3_time = nansum(tempvec2(1:3000))/sum(~isnan(tempvec2(1:3000)));
                rep_prop15_time = nansum(tempvec2(500:1500))/sum(~isnan(tempvec2(500:1500)));
            end
        else
            continue
        end
        
        
        %%
        %fixation duration by ordinal fixation number
        fixdurs = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        
        [recurrence_rate,recurrence_map,corm,laminarity,laminarity_len,...
            forward_trace,forward_trace_len,reverse_trace,reverse_trace_len] = ...
            calculate_auto_recurrence(nov_fix,N_fix,dist_threshold);
        
        img_id = img_id+1;
        all_img_types(img_id) = rep_trialtype;
        all_recurrence_rates(img_id) = recurrence_rate;
        all_laminarity(img_id) = laminarity;
        all_forward_trace(img_id) = forward_trace;
        all_reverse_trace(img_id) =reverse_trace;
        all_corm(img_id) = corm;
        
        
        all_ROI_time(img_id) = rep_propall_time;
        all_ROI_515(img_id) = rep_prop15_time;
        
        if sum(tempvec2) == 0
            all_ROI_first(img_id) = 0;
        else
            all_ROI_first(img_id) = min(find(tempvec2));
        end
        mean_fix_durs(img_id) = mean(fixdurs(4:end));
        
    end
end

%%
figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_forward_trace(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Fwd Retrace Rate')
    ylabel('% of Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(all_reverse_trace(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(all_reverse_trace(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Reverse Retrace Rate')
    ylabel('% of Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_forward_trace(all_img_types == type)+all_reverse_trace(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)'+all_reverse_trace(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('All Retrace Rate')
    ylabel('% of Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end

figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_recurrence_rates(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(all_recurrence_rates(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Recurrence Rate')
    ylabel('% of Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(mean_fix_durs(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(mean_fix_durs(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_laminarity(all_img_types == type),all_ROI_time(all_img_types == type),'k.')
    axis square
    rho =  corr(all_laminarity(all_img_types == type)',all_ROI_time(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end
%%
figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_forward_trace(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Fwd Retrace Rate')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(all_reverse_trace(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(all_reverse_trace(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Reverse Retrace Rate')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_forward_trace(all_img_types == type)+all_reverse_trace(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)'+all_reverse_trace(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('All Retrace Rate')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end

figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_recurrence_rates(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(all_recurrence_rates(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Recurrence Rate')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(mean_fix_durs(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(mean_fix_durs(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_laminarity(all_img_types == type),all_ROI_515(all_img_types == type),'k.')
    axis square
    rho =  corr(all_laminarity(all_img_types == type)',all_ROI_515(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of 0.5-1.5 Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end
%%
figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_forward_trace(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Fwd Retrace Rate')
    ylabel('% of Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(all_reverse_trace(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(all_reverse_trace(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Reverse Retrace Rate')
    ylabel('% of Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_forward_trace(all_img_types == type)+all_reverse_trace(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(all_forward_trace(all_img_types == type)'+all_reverse_trace(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('All Retrace Rate')
    ylabel('% of Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end

figure

for type = 2:4
    subplot(3,3,type-1)
    plot(all_recurrence_rates(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(all_recurrence_rates(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Recurrence Rate')
    ylabel('% of Time in ROI')
    
    subplot(3,3,type-1+3)
    plot(mean_fix_durs(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(mean_fix_durs(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of Time in ROI')
    
    
    subplot(3,3,type-1+6)
    plot(all_laminarity(all_img_types == type),all_ROI_first(all_img_types == type),'k.')
    axis square
    rho =  corr(all_laminarity(all_img_types == type)',all_ROI_first(all_img_types == type)','row','pairwise','type','Spearman');...
        title(['\rho = ' num2str(rho,2)])
    xlabel('Avg Fixation Duration (ms)')
    ylabel('% of Time in ROI')
    
    subtitle('Familiar/Replaced/Moved')
end
%%
figure
for type = 2:4
    group_durs = mean_fix_durs(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted into Good vs Bad by Fixation Duration')
%%

figure
for type = 2:4
    group_durs = all_recurrence_rates(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by Recurrence Rate')
%%
figure
for type = 2:4
    group_durs = all_laminarity(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by Laminarity Rate')
%%
figure
for type = 2:4
    group_durs = all_forward_trace(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by Forward Retrace Rate')
%%
figure
for type = 2:4
    group_durs = all_reverse_trace(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by Reverse Retrace Rate')
%%
figure
for type = 2:4
    group_durs = all_reverse_trace(all_img_types == type)+all_forward_trace(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by Whole Retrace Rate')
%%
figure
for type = 2:4
    group_durs = all_corm(all_img_types == type);
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    pct33 = prctile(group_durs,33);
    pct67 = prctile(group_durs,66);
    
    bad_ROI = group_ROI(group_durs <= pct33);
    bad_ROI15 = group_ROI15(group_durs <= pct33);
    
    good_ROI = group_ROI(group_durs >= pct67);
    good_ROI15 = group_ROI15(group_durs >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('Sorted by CORM')
%%
for type = 2:4
    group_vals = [all_recurrence_rates(all_img_types == type)' all_corm(all_img_types == type)'...
        all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...
        all_reverse_trace(all_img_types == type)' mean_fix_durs(all_img_types == type)'];
    group_ROI = all_ROI_time(all_img_types == type);
    group_ROI15 = all_ROI_515(all_img_types == type);
    
    T = kmeans(group_vals,3);
    
    
    %         Y = pdist(group_vals);
    %     Z = linkage(Y);
    %     dendrogram(Z)
    %
    %     T = cluster(Z,'Depth',9)
    
    [U,S,V] = pca(group_vals,1);
    pct33 = prctile(U,33);
    pct67 = prctile(U,67);
    
    bad_ROI = group_ROI(U <= pct33);
    bad_ROI15 = group_ROI15(U <= pct33);
    
    good_ROI = group_ROI(U >= pct67);
    good_ROI15 = group_ROI15(U >= pct67);
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    [~,p15] = kstest2(good_ROI15,bad_ROI15);
    
    subplot(2,2,type-1)
    hold on
    bar([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]])
    errorb([[mean(bad_ROI) mean(good_ROI)];[mean(bad_ROI15) mean(good_ROI15)]],...
        [[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))];...
        [std(bad_ROI15)/sqrt(length(bad_ROI15)) std(good_ROI15)/sqrt(length(good_ROI15))]])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    legend('Whole Trial','0.5-1.5 s')
    ylabel('Time in ROI')
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2) ', p_{15} = '  num2str(p15,2)])
    end
end
subtitle('PCA sorted')
%%