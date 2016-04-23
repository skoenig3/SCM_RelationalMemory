function combinedPupilROI_SCMRM(data_dir,cortex_files,ROI_dir)


all_nov_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rep_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rpl_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_mvd_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples

nov_vals = [];
rep_vals = [];
rpl_vals = [];
mvd_vals = [];
for file = 1:length(cortex_files)
    %---Import Data--%
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(10) '-fixation.mat']);
    load([ROI_dir itemfile(1:end-4) '_ROIs.mat']);
    setnum = str2double(itemfile(6:7));
    
    %---Preallocate space for data structure---%
    allnumFixationsInROI = cell(1,4); %number of fixations in ROI by type
    num10FixationsInROI = cell(1,4); %number of fixations in ROI by type for 1st 10 fixations
    allTime_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type
    Time3_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type in 1st 3 secs
    Time15_ROI_timeWindow =  cell(1,4);%amount of time in ROI  by type in 0.5-1.5 secs window
    area = cell(1,4); %area of ROI by type
    
    
    nov_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    rep_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    rpl_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    mvd_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    
    
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
        
        nov_x = fixationstats{2*img-1}.XY(1,:);
        nov_y = fixationstats{2*img-1}.XY(2,:);
        nov_pupil = pupildata{2*img-1}(round(nov_img_on/5):round(nov_img_off/5));
        
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
        
        rep_x = fixationstats{2*img}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{2*img}.XY(2,rep_img_on:rep_img_off);
        rep_pupil = pupildata{2*img}(round(rep_img_on/5):round(rep_img_off/5));
        
        
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
        
        %---Find When the monkey Blinked,when monkey looked away---%
        %pupil values at 0 diameter
        
        [nov_blink_ind,nov_nan_ind,nov_time_out] = findbad_eye_ind(nov_pupil,nov_x,1000);
        [rep_blink_ind,rep_nan_ind,rep_time_out] = findbad_eye_ind(rep_pupil,rep_x,1000);
        
        %---Get Pupil Data over time---%
        nov_pupil = pupildata{2*img-1}((round(nov_img_on/5)-100):round(nov_img_off/5)); %grab 500 ms before img on too
        rep_pupil = pupildata{2*img}((round(rep_img_on/5)-100):round(rep_img_off/5)); %grab 500 ms before img on too
        
        %remove blinks from pupil data
        if ~isempty(nov_blink_ind)
            for b = 1:size(nov_blink_ind,1)
                ind = nov_blink_ind(b,:);
                ind(ind == 0) = [];
                nov_pupil(ind+100) = NaN;
            end
        end
        if ~isempty(rep_blink_ind)
            for b = 1:size(rep_blink_ind,1)
                ind = rep_blink_ind(b,:);
                ind(ind == 0) = [];
                rep_pupil(ind+100) = NaN;
            end
        end
        
        %remove pupil data when not looking at picture
        if ~isempty(nov_nan_ind)
            for b = 1:size(nov_nan_ind,1)
                ind = nov_nan_ind(b,:);
                ind(ind == 0) = [];
                if ind(1) < 5
                    ind(1) = 5;
                end
                nov_pupil(round(ind(1)/5):round(ind(end)/5)) = NaN;
            end
        end
        if ~isempty(rep_nan_ind)
            for b = 1:size(rep_nan_ind,1)
                ind = rep_nan_ind(b,:);
                ind(ind == 0) = [];
                if ind(1) < 5
                    ind(1) = 5;
                end
                rep_pupil(round(ind(1)/5):round(ind(end)/5)) = NaN;
            end
        end
        
        %remove data after looked away too much
        if ~isnan(nov_time_out)
            nov_pupil(round(nov_time_out/5):end) = NaN;
        end
        if ~isnan(rep_time_out)
            rep_pupil(round(rep_time_out/5):end) = NaN;
        end
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)))+2;
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)))+2;
        
        %only look at the 1st 7 seconds the images is up
        nov_pupil = nov_pupil(1:1500);
        rep_pupil = rep_pupil(1:1500);
        
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
            if ~isempty(ROIfix1)
                area_ROI1 = abs(ROI1(1)-ROI1(2))*abs(ROI1(3)-ROI1(4)); %area of ROI for novel image
                area_ROI2 = abs(ROI2(1)-ROI2(2))*abs(ROI2(3)-ROI2(4)); %for replaced, moved or familiar
                
                area{1} = [area{1} area_ROI1]; %area of novel
                area{rep_trialtype} = [area{rep_trialtype} area_ROI2];%area of replaced, moved or familiar
            end
            
            
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
            
                        
            nov_trial_by_trial_pupil(img,:) = nov_pupil;
            nov_vals = [nov_vals [nanmean(nov_pupil(300:600)); nov_prop3_time]];
            if rep_trialtype == 2 %repeat
                rep_trial_by_trial_pupil(img,:) = rep_pupil;
                rep_vals = [rep_vals [nanmean(rep_pupil(300:600)); rep_prop3_time]];
            elseif rep_trialtype == 3 %replaced
                rpl_trial_by_trial_pupil(img,:) = rep_pupil;
                rpl_vals = [rpl_vals [nanmean(rep_pupil(300:600)); rep_prop3_time]];
            elseif rep_trialtype == 4 %moved
                mvd_trial_by_trial_pupil(img,:) = rep_pupil;
                mvd_vals = [mvd_vals [nanmean(rep_pupil(300:600)); rep_prop3_time]];
            else
                error('unknown 2nd presentation type')
            end
        end        
    end
    
    all_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    all_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    all_rpl_pupil(file,:) = nanmean(rpl_trial_by_trial_pupil);
    all_mvd_pupil(file,:) = nanmean(mvd_trial_by_trial_pupil);
end
%%
%---Plot Pupil Diameter Over Time---%
t = 1:5:7500;
figure
hold on
dofill(t,all_nov_pupil/1000,'blue',1,60)
dofill(t,all_rep_pupil/1000,'red',1,60)
dofill(t,all_rpl_pupil/1000,'magenta',1,60)
dofill(t,all_mvd_pupil/1000,'green',1,60)
yl = ylim;
plot([500 500],[0.75 1.1],'--k')
hold off
ylim(yl)
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
legend('Novel','Repeat','Relaced','Moved')
title([cortex_files{1}(1:2) ' : Normalized Pupil Horizontal Diameter'])
%%
figure
subplot(2,2,1)
plot(nov_vals(1,:),100*nov_vals(2,:),'.')
xlabel('Pupil Diameter')
ylabel('% of 0.5-1.5 sec in ROI')
title('Novel')
ylim([-5 105])

subplot(2,2,2)
plot(rep_vals(1,:),100*rep_vals(2,:),'.')
xlabel('Pupil Diameter')
ylabel('% of 0.5-1.5 sec in ROI')
title('Repeats')
ylim([-5 105])


subplot(2,2,3)
plot(rpl_vals(1,:),100*rpl_vals(2,:),'.')
xlabel('Pupil Diameter')
ylabel('% of 0.5-1.5 sec in ROI')
title('Replaced')
ylim([-5 105])


subplot(2,2,4)
plot(mvd_vals(1,:),100*mvd_vals(2,:),'.')
xlabel('Pupil Diameter')
ylabel('% of 0.5-1.5 sec in ROI')
title('Moved')
ylim([-5 105])
subtitle([cortex_files{1}(1:2) ': Relationship of Behavior and Pupil'])
disp('done')
