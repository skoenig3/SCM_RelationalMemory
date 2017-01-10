function SCMRM_ROIanalysis(cortexfile,itemfile,figure_dir,data_dir,ROI_dir,img_dir)
% written by Seth Koenig 9/1/15. Modified from SCMRM_analysis.m and VR_SCM_analysis.m
% code process main analysis for the SCMRM task and determines the number
% of fixations and the amount of time spent in the ROI. ROIs are determined
% with computer assistance using scmgui_RM.m.

ploteyetraces = true;%set to false if don't want to plot eye traces on individual images

%---Import Data--%
load([data_dir cortexfile(1:8) '_' cortexfile(10) '-fixation.mat']);
load([ROI_dir itemfile(1:end-4) '_ROIs.mat']);
figure_dir2 = [figure_dir itemfile(1:end-4) '\']; %where to save individual eye traces and other set plots
mkdir(figure_dir2)
setnum = str2double(itemfile(6:7));

%---Preallocate space for data structure---%
allnumFixationsInROI = cell(1,4); %number of fixations in ROI by type
num10FixationsInROI = cell(1,4); %number of fixations in ROI by type for 1st 10 fixations
allTime_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type
Time3_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type in 1st 3 secs
Time15_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type in 0.5-1.5 secs window
area = cell(1,4); %area of ROI by type
image_numbers = cell(1,4); %image number to help keep track

%time in ROI over time
all_time = zeros(5,7000);%novel, familiar, replaced, moved,new_moved
all_time_points = zeros(5,7000);

setname = ['SCM' itemfile(6:7)];

imageX = 800; %horiztonal image size
imageY = 600; %vertical image size
for img = 1:36
    
    if img == 18 && setnum == 30%image accidentaly was shown in a previous set
       continue 
    end
    
    %---get image names---%
    nov_name = image_names{1,img};
    rep_name = image_names{2,img};
    rep_trialtype = trialtype(2,img);
    
    %---Grab Important Vairables---%
    %for novel image
    nov_allval = per(2*img-1).allval;
    nov_alltim = per(2*img-1).alltim;
    nov_img_on = nov_alltim(nov_allval == 23)-nov_alltim(nov_allval == 100);%image on relative to eye data start
    nov_img_off = nov_alltim(nov_allval == 24)-nov_alltim(nov_allval == 100);%image on relative to eye data start
   
    if nov_img_off-nov_img_on > 20000 %EOG overlfow. Probably ony RED
       continue 
    end
    
    nov_x = fixationstats{2*img-1}.XY(1,:);
    nov_y = fixationstats{2*img-1}.XY(2,:);
    
    nov_fix = fixationstats{2*img-1}.fixations;
    nov_fixtimes = fixationstats{2*img-1}.fixationtimes;
    
    pre_img_fix = find(nov_fixtimes(1,:) <= nov_img_on);
    post_img_fix = find(nov_fixtimes(2,:) > nov_img_off);
    
    if ~isempty(post_img_fix)
        nov_x(nov_fixtimes(1,post_img_fix(1)):end) = [];
        nov_y(nov_fixtimes(1,post_img_fix(1)):end) = [];
    end
    if ~isempty(pre_img_fix)
        nov_x(1:nov_fixtimes(2,pre_img_fix(end))) = [];
        nov_y(1:nov_fixtimes(2,pre_img_fix(end))) = [];
    end
    
    nov_fix(:,post_img_fix) = [];
    nov_fix(:,pre_img_fix) = [];
    nov_fixtimes(:,post_img_fix) = [];
    nov_fixtimes(:,pre_img_fix ) = [];
    
    nov_fixtimes = nov_fixtimes-nov_img_on; %center to image start
    
    %for repeat image
    rep_allval = per(2*img).allval;
    rep_alltim = per(2*img).alltim;
    rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
    rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
    
    if rep_img_off-rep_img_on > 20000 %EOG overlfow
        continue
    end
    
    rep_x = fixationstats{2*img}.XY(1,:);
    rep_y = fixationstats{2*img}.XY(2,:);
    
    rep_fix = fixationstats{2*img}.fixations;
    rep_fixtimes = fixationstats{2*img}.fixationtimes;
    
    pre_img_fix = find(rep_fixtimes(1,:) <= rep_img_on);
    post_img_fix = find(rep_fixtimes(2,:) > rep_img_off);
    
    if ~isempty(post_img_fix)
        rep_x(rep_fixtimes(1,post_img_fix(1)):end) = [];
        rep_y(rep_fixtimes(1,post_img_fix(1)):end) = [];
    end
    if ~isempty(pre_img_fix)
        rep_x(1:rep_fixtimes(2,pre_img_fix(end))) = [];
        rep_y(1:rep_fixtimes(2,pre_img_fix(end))) = [];
    end
    
    rep_fix(:,post_img_fix) = [];
    rep_fix(:,pre_img_fix) = [];
    rep_fixtimes(:,post_img_fix) = [];
    rep_fixtimes(:,pre_img_fix ) = [];
    
    rep_fixtimes = rep_fixtimes-rep_img_on;%center at images start
    
    %---Remove data outside image and afte 7 seconds---%
    nov_time_out = zeros(1,7000); %time outside image or crosshair fixation
    nov_outside = find(isnan(nov_x));
    nov_outside(nov_outside > 7000) = [];
    nov_time_out(nov_outside) = 1;
    late = find(nov_fixtimes(1,:) > 7000); %ignore anything after 7000 ms
    nov_fixtimes(:,late) = [];
    nov_fix(:,late) = [];
    nov_fixtimes(nov_fixtimes > 7000) = 7000; %and cap end times
    
    rep_time_out = zeros(1,7000); %time outside image or crosshair fixation
    rep_outside = find(isnan(rep_x));
    rep_outside(rep_outside > 7000) = [];
    rep_time_out(rep_outside) = 1;
    late = find(rep_fixtimes(1,:) > 7000); %ignore anything after 7000 ms
    rep_fixtimes(:,late) = [];
    rep_fix(:,late) = [];
    rep_fixtimes(rep_fixtimes > 7000) = 7000; %and cap end times
    
    
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
        
        area{1} = [area{1} area_ROI1]; %area of novel
        area{rep_trialtype} = [area{rep_trialtype} area_ROI2];%area of replaced, moved or familiar
        
        %store image numbers so data is easier to track
        image_numbers{1} = [image_numbers{1} img];
        image_numbers{rep_trialtype} = [image_numbers{rep_trialtype} img];

        %---Calculate Average Time in ROI---%
        tempvec1=zeros(1,7000);%novel image
        for j=1:length(ROIfix1)
            tempvec1(nov_fixtimes(1,ROIfix1(j)):nov_fixtimes(2,ROIfix1(j)))=1;
        end
        tempvec1(find(nov_time_out)) = NaN;
        
        tempvec2=zeros(1,7000);%manipulated image new ROI
        for j=1:length(ROIfix2)
            tempvec2(rep_fixtimes(1,ROIfix2(j)):rep_fixtimes(2,ROIfix2(j)))=1;
        end
        tempvec2(find(rep_time_out)) = NaN;
        
        tempvec3=zeros(1,7000);%manipulated imaged old moved ROI
        for j=1:length(ROIfix2_oldROI)
            tempvec3(rep_fixtimes(1,ROIfix2_oldROI(j)):rep_fixtimes(2,ROIfix2_oldROI(j)))=1;
        end
        tempvec3(find(rep_time_out)) = NaN;
        
        
        %---Calculate Stats for the Novel image---%
        nov_propall_fixations = length(ROIfix1)/length(nov_fix);
        nov_prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
        nov_propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
        nov_prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
        nov_prop15_time = nansum(tempvec1(500:1500))/sum(~isnan(tempvec1(500:1500)));
        
        allnumFixationsInROI{1} =[allnumFixationsInROI{1} nov_propall_fixations];
        num10FixationsInROI{1} = [num10FixationsInROI{1} nov_prop10_fixaitons];
        allTime_ROI_timeWindow{1} = [allTime_ROI_timeWindow{1} nov_propall_time];
        Time3_ROI_timeWindow{1} = [Time3_ROI_timeWindow{1} nov_prop3_time];
        Time15_ROI_timeWindow{1}= [ Time15_ROI_timeWindow{1} nov_prop15_time];
        
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
            
            allnumFixationsInROI{rep_trialtype} =[allnumFixationsInROI{rep_trialtype} rep_propall_fixations];
            num10FixationsInROI{rep_trialtype} = [num10FixationsInROI{rep_trialtype} rep_prop10_fixaitons];
            allTime_ROI_timeWindow{rep_trialtype} = [allTime_ROI_timeWindow{rep_trialtype} rep_propall_time];
            Time3_ROI_timeWindow{rep_trialtype} = [Time3_ROI_timeWindow{rep_trialtype} rep_prop3_time];
            Time15_ROI_timeWindow{rep_trialtype} = [Time15_ROI_timeWindow{rep_trialtype} rep_prop15_time];
            
        else %for replaced or repeat
            rep_propall_fixations = length(ROIfix2)/length(rep_fix);
            rep_prop10_fixaitons = length(ROIfix2(ROIfix2 <= 10))/10;
            rep_propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
            rep_prop3_time = nansum(tempvec2(1:3000))/sum(~isnan(tempvec2(1:3000)));
            rep_prop15_time = nansum(tempvec2(500:1500))/sum(~isnan(tempvec2(500:1500)));
            
            allnumFixationsInROI{rep_trialtype} =[allnumFixationsInROI{rep_trialtype} rep_propall_fixations];
            num10FixationsInROI{rep_trialtype} = [num10FixationsInROI{rep_trialtype} rep_prop10_fixaitons];
            allTime_ROI_timeWindow{rep_trialtype} = [allTime_ROI_timeWindow{rep_trialtype} rep_propall_time];
            Time3_ROI_timeWindow{rep_trialtype} = [Time3_ROI_timeWindow{rep_trialtype} rep_prop3_time];
            Time15_ROI_timeWindow{rep_trialtype} = [Time15_ROI_timeWindow{rep_trialtype} rep_prop15_time];
        end
        
        
        %---Combined Time in ROI over Time across images by type---%
        all_time_points(1,~isnan(tempvec1)) = all_time_points(1,~isnan(tempvec1))+1;
        tempvec1(isnan(tempvec1)) = 0;
        all_time(1,:) = all_time(1,:) + tempvec1;
        
        all_time_points(rep_trialtype,~isnan(tempvec2)) = ...
            all_time_points(rep_trialtype,~isnan(tempvec2))+1;
        tempvec2(isnan(tempvec2)) = 0;
        all_time(rep_trialtype,:) =  all_time(rep_trialtype,:)+tempvec2;
        
        if rep_trialtype == 4 %for moved only
            all_time_points(5,~isnan(tempvec3)) = ...
                all_time_points(5,~isnan(tempvec3))+1;
            tempvec3(isnan(tempvec3)) = 0;
            all_time(5,:) =  all_time(5,:)+tempvec3+tempvec2;
        end
    end
    
    %---Visual Verification and Plt Eye Traces by Image---%
    if ploteyetraces
        figure
        subplot(2,2,1)
        hold on
        image(flipud(imread([img_dir setname '\' nov_name])));
        plot(nov_x,nov_y);
        plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
            [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
        plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
            [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
        hold off
        box off
        xlim([0 imageX])
        ylim([0 imageY])
        axis off
        axis equal
        title(nov_name)
        
        subplot(2,2,3)
        hold on
        image(flipud(imread([img_dir setname '\' rep_name])));
        plot(rep_x,rep_y);
        plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
            [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
        plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
            [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
        hold off
        box off
        xlim([0 imageX])
        ylim([0 imageY])
        axis off
        axis equal
        title(rep_name)
        
        if isempty(ROIfix1)%didnt look at original region, moving on to next image
            subplot(2,2,2)
            title('No fixations in Novel ROI')
            save_and_close_fig(figure_dir2,[cortexfile(1:2) '_' num2str(img)]);
            continue
        end
        
        s(1) = subplot(2,2,2);
        bar([nov_propall_fixations nov_propall_time nov_prop10_fixaitons nov_prop3_time])
        set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
        ylabel('Proportion in ROI')
        title(['First Fixation in ROI is # ' num2str(ROIfix1(1))])
        yl1 = ylim;
        
        s(2) = subplot(2,2,4);
        bar([rep_propall_fixations rep_propall_time rep_prop10_fixaitons rep_prop3_time])
        set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
        ylabel('Proportion in ROI')
        if ~isempty(ROIfix2)
            title(['First Fixation in ROI is # ' num2str(ROIfix2(1))])
        end
        yl2 = ylim;
        
        ymax = max(yl2(2),yl1(2));
        subplot(2,2,2)
        ylim([0 ymax])
        subplot(2,2,4)
        ylim([0 ymax])
        
        
        
        subtitle(cortexfile(1:2))
        save_and_close_fig(figure_dir2,[cortexfile(1:2) '_' num2str(img)]);
    end
    
end


%---Plot Result Averaged across the Image Set---%
%get average time over time
all_time_points(all_time_points == 0) = 1; %for division purposes
all_time = all_time./all_time_points; %normalize

sess_num_fix =NaN(1,4); %proportion of all fixations in ROI
sess_num10fix = NaN(1,4); %proportion of 1st 10 fixations in ROI
sess_allTime = NaN(1,4); %proportion of all time in ROI
sess_Time3 = NaN(1,4); %proportion of 1st 3 seconds in ROI
sess_Time15 =  NaN(1,4);%proportion of 0.5-1.5 seconds in ROI
for tp = 1:4 %for each trial type
    sess_num_fix(tp) = nanmean(allnumFixationsInROI{tp});
    sess_num10fix(tp) = nanmean(num10FixationsInROI{tp});
    sess_allTime(tp) = nanmean(allTime_ROI_timeWindow{tp});
    sess_Time3(tp) = nanmean(Time3_ROI_timeWindow{tp});
    sess_Time15(tp)= nanmean(Time15_ROI_timeWindow{tp});
end

figure
subplot(2,3,1)
bar(sess_num10fix);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of 1st 10 Fixations in ROI')

subplot(2,3,2)
bar(sess_Time3);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of Time in 1st 3 seconds in ROI')

subplot(2,3,3)
bar(sess_num_fix);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Fixations in ROI')

subplot(2,3,4)
bar(sess_allTime);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Time in ROI')

subplot(2,3,5)
bar(sess_Time15);
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of Time in 0.5-1.5 sec window in ROI')

clrs = 'brmgg';
line = {'-','-','-','-','--'};
subplot(2,3,6)
hold on
for r = 1:size(all_time,1)
    plot(100*filtfilt(1/200*ones(1,200),1,all_time(r,:)),...
        [clrs(r) line{r}]);
end
hold off
xlabel('Time (ms)')
ylabel('% of time in ROI')
legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
title('% of Time in ROI over Time')

save_and_close_fig(figure_dir2,[cortexfile(1:2) '_' setname]);


save([data_dir cortexfile(1:8) '_' cortexfile(10) '-ROIdata.mat'],...
    'allnumFixationsInROI','num10FixationsInROI','allTime_ROI_timeWindow',...
    'Time3_ROI_timeWindow','Time15_ROI_timeWindow','area','all_time',...
    'image_numbers','itemfile');