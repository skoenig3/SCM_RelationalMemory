function combinedEyeMovementStats_SCMRM(data_dir,cortex_files)
%function combines viewing behavior across multiple sessions
% written by Seth Koenig
% 1) calculate fixation durations
% 2) calculate saccade amplitudes
% 3) calculate pupil diameter 

first_fix = cell(1,length(cortex_files));
first_fix_location = cell(2,length(cortex_files));
all_blinks = cell(1,2);


all_nov_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rep_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_rpl_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_mvd_pupil = NaN(length(cortex_files),1500);%7 seconds sampled at 200 Hz is 1400 samples
all_novel_fix_dur = NaN(length(cortex_files),50); %fixation durations
all_repeat_fix_dur = NaN(length(cortex_files),50); %fixation durations
all_novel_sac_amp = NaN(length(cortex_files),50); %saccade amplitudes
for file = 1:length(cortex_files)
    load([data_dir cortex_files{file}(1:8) '_' cortex_files{file}(end) '-fixation.mat'])
    setnum = str2double(itemfile(6:7));
    
    
    nov_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    rep_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    rpl_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    mvd_trial_by_trial_pupil = NaN(36,1500); %each trials pupil data
    nov_trial_by_trial_fixdurs = NaN(36,50); %each trials fixation durations
    rep_trial_by_trial_fixdurs = NaN(36,50); %each trials fixation durations
    nov_trial_by_trial_sacamps = NaN(36,50); %each trials saccade amplitudes
    rep_trial_by_trial_sacamps = NaN(36,50); %each trials fixation amplitudes
    
    
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
        
        rep_x = fixationstats{2*img}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{2*img}.XY(2,rep_img_on:rep_img_off);
        rep_pupil = pupildata{2*img}(round(rep_img_on/5):round(rep_img_off/5));
        
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
                nov_pupil(round(ind(1)/5):round(ind(end)/5)) = NaN;
            end
        end
        if ~isempty(rep_nan_ind)
            for b = 1:size(rep_nan_ind,1)
                ind = rep_nan_ind(b,:);
                ind(ind == 0) = [];
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
        
        nov_trial_by_trial_pupil(img,:) = nov_pupil;
        if trialtype(2,img) == 2 %repeat
            rep_trial_by_trial_pupil(img,:) = rep_pupil;
        elseif trialtype(2,img) == 3 %replaced
            rpl_trial_by_trial_pupil(img,:) = rep_pupil;
        elseif trialtype(2,img) == 4 %moved
            mvd_trial_by_trial_pupil(img,:) = rep_pupil;
        else
            error('unknown 2nd presentation type')
        end
        
        %---Calculate Fixation Durations and Saccade Amplitudes---%         
        if ~isnan(nov_time_out)
            post_attention_fix = find(nov_fixtimes(1,:) > nov_time_out);
            nov_fixtimes(:,post_attention_fix) = [];
        end
        
        if ~isnan(nov_time_out)
            post_attention_fix = find(rep_fixtimes(1,:) > rep_time_out);
            rep_fixtimes(:,post_attention_fix) = [];
        end
 
        %fixation duration by ordinal fixation number
        nov_trial_by_trial_fixdurs(img,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        rep_trial_by_trial_fixdurs(img,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
        nov_x = fixationstats{2*img-1}.XY(1,:);
        nov_y = fixationstats{2*img-1}.XY(2,:);
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            nov_trial_by_trial_sacamps(img,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = fixationstats{2*img}.XY(1,:);
        rep_y = fixationstats{2*img}.XY(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            rep_trial_by_trial_sacamps(img,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    all_nov_pupil(file,:) = nanmean(nov_trial_by_trial_pupil);
    all_rep_pupil(file,:) = nanmean(rep_trial_by_trial_pupil);
    all_rpl_pupil(file,:) = nanmean(rpl_trial_by_trial_pupil);
    all_mvd_pupil(file,:) = nanmean(mvd_trial_by_trial_pupil);
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_fixdurs'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    all_novel_fix_dur(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_fixdurs(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_fixdurs'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    all_repeat_fix_dur(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_fixdurs(:,1:rep_median_num_fix));
    
    
    nov_median_num_fix = sum(~isnan(nov_trial_by_trial_sacamps'));
    nov_median_num_fix(nov_median_num_fix == 0) = [];
    nov_median_num_fix = round(median(nov_median_num_fix));
    all_novel_sac_amp(file,1:nov_median_num_fix) = ...
        nanmean(nov_trial_by_trial_sacamps(:,1:nov_median_num_fix));
    
    rep_median_num_fix = sum(~isnan(rep_trial_by_trial_sacamps'));
    rep_median_num_fix(rep_median_num_fix == 0) = [];
    rep_median_num_fix = round(median(rep_median_num_fix));
    all_repeat_sac_amp(file,1:rep_median_num_fix) = ...
        nanmean(rep_trial_by_trial_sacamps(:,1:rep_median_num_fix));
% 
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
%---Plot Fixation Durations By Fixation Number---%
nov_median_num_fix = ceil(median(sum(~isnan(all_novel_fix_dur'))));
rep_median_num_fix = ceil(median(sum(~isnan(all_repeat_fix_dur'))));
all_novel_fix_dur = all_novel_fix_dur(:,1:nov_median_num_fix);
all_repeat_fix_dur = all_repeat_fix_dur(:,1:rep_median_num_fix);

figure
hold on
errorbar(nanmean(all_novel_fix_dur),nanstd(all_novel_fix_dur)...
    ./sqrt(sum(~isnan(all_novel_fix_dur))),'b')
errorbar(nanmean(all_repeat_fix_dur),nanstd(all_repeat_fix_dur)...
    ./sqrt(sum(~isnan(all_repeat_fix_dur))),'r')
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))
%%
%---Plot Saccade Ampltiudes By Saccades Number---%
nov_median_num_sac = ceil(median(sum(~isnan(all_novel_sac_amp'))));
rep_median_num_sac = ceil(median(sum(~isnan(all_repeat_sac_amp'))));
all_novel_sac_amp = all_novel_sac_amp(:,1:nov_median_num_sac)/24;
all_repeat_sac_amp = all_repeat_sac_amp(:,1:rep_median_num_sac)/24;

figure
hold on
errorbar(nanmean(all_novel_sac_amp),nanstd(all_novel_sac_amp)...
    ./sqrt(sum(~isnan(all_novel_sac_amp))),'b')
errorbar(nanmean(all_repeat_sac_amp),nanstd(all_repeat_sac_amp)...
    ./sqrt(sum(~isnan(all_repeat_sac_amp))),'r')
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Novel','Repeat')
title(cortex_files{1}(1:2))