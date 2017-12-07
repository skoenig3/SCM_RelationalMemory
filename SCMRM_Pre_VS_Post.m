%compare pre-post lesion
%written by Seth Konig 6/7/16

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';%location of processed eye data
ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\'; %locations of regions of interest
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Figures\'; %where to put plots
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\'; %where all the image sets are

%---Vivian---%

% pre_files = {'PW150901.2','PW150902.2','PW150903.2','PW150904.2',...
%     'PW150908.2','PW150910.2','PW150911.2',...
%     'PW150914.2','PW150916.2','PW150917.2','PW150918.2',...
%     'PW150921.2','PW150922.2','PW150923.2','PW150924.2'};
%
% post_files = {'PW160401.2','PW160405.2','PW160406.2','PW160407.2','PW160408.2',...
%     'PW160411.2','PW160412.2','PW160413.2','PW160414.2','PW160415.2',...
%     'PW160418.2','PW160419.2','PW160420.2','PW160421.2','PW160422.2'};

%---Red---%
% pre_files =      {'RR151027.2','RR151026.2','RR151023.2','RR151022.2',...
%     'RR151021.2','RR151020.2','RR151019.2','RR151016.2',...
%     'RR151015.3','RR151014.2','RR151013.2','RR151009.2',...
%     'RR151008.2','RR151007.2','RR151005.2','RR151002.2',...
%     'RR151001.2'};
% post_files = {'RR160425.2','RR160426.2','RR160701.2','RR160706.2','RR160707.2',...
%     'RR160708.2','RR160714.2','RR160715.2','RR160718.2','RR160721.2',...
%     'RR160725.2','RR160726.2','RR160727.2','RR160728.2','RR160729.2',...
%     'RR160801.2','RR160803.2','RR160804.2','RR160805.2'};
%

%---Tobii---%
pre_files = {'TO150901.2','TO150903.2','TO150904.3',...
    'TO150908.2','TO150909.2','TO150910.2','TO150911.2',...
    'TO150915.2','TO150916.2','TO150918.2',...
    'TO151001.2','TO151002.2','TO151006.2','TO151007.2'};

post_files = {'TO170620.2','TO170621.2','TO170622.2','TO170623.2',...
    'TO170626.2','TO170627.2','TO170628.2','TO170629.2',...
    'TO170630.2','TO170705.2','TO170706.2','TO170707.2',...
    'TO170710.2','TO170711.2','TO170712.2'};


%---Get Pre-Data---%
pre_num_sets = length(pre_files);
pre_sets_numFix = NaN(length(pre_files),4);
pre_sets_10numFix = NaN(length(pre_files),4);
pre_sets_Time = NaN(length(pre_files),4);
pre_sets_Time3 = NaN(length(pre_files),4);
pre_sets_Time15 = NaN(length(pre_files),4);
pre_time_sets = zeros(5,7000);

%Eye data
pre_pupil = cell(1,4);
pre_sacamp = cell(1,4);
pre_fixdur = cell(1,4);

for sets = 1:length(pre_files)
    %---ROI anlalysis---%
    load([data_dir pre_files{sets}(1:8) '_' pre_files{sets}(10) '-ROIdata.mat'])
    
    for type = 1:4
        num_points(sets,type)= sum(~isnan(allnumFixationsInROI{type}));
        pre_sets_numFix(sets,type) = 100*nanmean(allnumFixationsInROI{type});
        pre_sets_10numFix(sets,type) = 100*nanmean(num10FixationsInROI{type});
        pre_sets_Time(sets,type) = 100*nanmean(allTime_ROI_timeWindow{type});
        pre_sets_Time3(sets,type) = 100*nanmean(Time3_ROI_timeWindow{type});
        pre_sets_Time15(sets,type) = 100*nanmean(Time15_ROI_timeWindow{type});
    end
    pre_time_sets = pre_time_sets +all_time;
    
    %---Eye Movement Analysis---%
    load([data_dir pre_files{sets}(1:8) '_' pre_files{sets}(10) '-fixation.mat'])
    setnum = str2double(itemfile(6:7));
    
    trial_types = NaN(1,72);
    pupil_data = NaN(72,1500); %each trials pupil data
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
        
        pupil_data(2*img-1,:) = nov_pupil(1:1500);
        pupil_data(2*img,:) = rep_pupil(1:1500);
        
        trial_types(2*img-1)= 1;
        if trialtype(2,img) == 2 %repeat
            trial_types(2*img) = 2;
        elseif trialtype(2,img) == 3 %replaced
            trial_types(2*img) = 3;
        elseif trialtype(2,img) == 4 %moved
            trial_types(2*img) = 4;
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
        fixdurs(2*img-1,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        fixdurs(2*img,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
        nov_x = fixationstats{2*img-1}.XY(1,:);
        nov_y = fixationstats{2*img-1}.XY(2,:);
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            sacamps(2*img-1,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = fixationstats{2*img}.XY(1,:);
        rep_y = fixationstats{2*img}.XY(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            sacamps(2*img,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    %---Store data across all sets---%
    if size(sacamps,2) > 50
        sacamps = sacamps(:,1:50);
    end
    if size(fixdurs,2) > 50
        fixdurs = fixdurs(:,1:50);
    end
    for type = 1:4
        if sum(type == trial_types) > 1
            %pre_pupil{type} = [pre_pupil{type}; nanmean(pupil_data(trial_types == type,:))];
            pre_sacamp{type} = [pre_sacamp{type}; nanmean(sacamps(trial_types == type,:))];
            pre_fixdur{type} = [pre_fixdur{type}; nanmean(fixdurs(trial_types == type,:))];
        end
    end
end
pre_time_sets = pre_time_sets/pre_num_sets;

%---Get Post-Data---%
post_num_sets = length(post_files);
post_sets_numFix = NaN(length(post_files),4);
post_sets_10numFix = NaN(length(post_files),4);
post_sets_Time = NaN(length(post_files),4);
post_sets_Time3 = NaN(length(post_files),4);
post_sets_Time15 = NaN(length(post_files),4);
post_time_sets = zeros(5,7000);

%Eye data
post_pupil = cell(1,4);
post_sacamp = cell(1,4);
post_fixdur = cell(1,4);

for sets = 1:length(post_files)
    %---ROI anlalysis---%
    load([data_dir post_files{sets}(1:8) '_' post_files{sets}(10) '-ROIdata.mat'])
    
    for type = 1:4
        num_points(sets,type)= sum(~isnan(allnumFixationsInROI{type}));
        post_sets_numFix(sets,type) = 100*nanmean(allnumFixationsInROI{type});
        post_sets_10numFix(sets,type) = 100*nanmean(num10FixationsInROI{type});
        post_sets_Time(sets,type) = 100*nanmean(allTime_ROI_timeWindow{type});
        post_sets_Time3(sets,type) = 100*nanmean(Time3_ROI_timeWindow{type});
        post_sets_Time15(sets,type) = 100*nanmean(Time15_ROI_timeWindow{type});
    end
    post_time_sets = post_time_sets +all_time;
    
    %---Eye Movement Analysis---%
    load([data_dir post_files{sets}(1:8) '_' post_files{sets}(10) '-fixation.mat'])
    setnum = str2double(itemfile(6:7));
    
    trial_types = NaN(1,72);
    pupil_data = NaN(72,1500); %each trials pupil data
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
        
        pupil_data(2*img-1,:) = nov_pupil(1:1500);
        pupil_data(2*img,:) = rep_pupil(1:1500);
        
        trial_types(2*img-1)= 1;
        if trialtype(2,img) == 2 %repeat
            trial_types(2*img) = 2;
        elseif trialtype(2,img) == 3 %replaced
            trial_types(2*img) = 3;
        elseif trialtype(2,img) == 4 %moved
            trial_types(2*img) = 4;
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
        fixdurs(2*img-1,1:length(nov_fixtimes)) = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        fixdurs(2*img,1:length(rep_fixtimes)) = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
        nov_x = fixationstats{2*img-1}.XY(1,:);
        nov_y = fixationstats{2*img-1}.XY(2,:);
        for sac = 1:size(nov_sactimes,2)
            sacx = nov_x(nov_sactimes(2,sac)) - nov_x(nov_sactimes(1,sac));
            sacy = nov_y(nov_sactimes(2,sac)) - nov_y(nov_sactimes(1,sac));
            sacamps(2*img-1,sac) = sqrt(sacx^2+sacy^2);
        end
        
        rep_x = fixationstats{2*img}.XY(1,:);
        rep_y = fixationstats{2*img}.XY(2,:);
        for sac = 1:size(rep_sactimes,2)
            sacx = rep_x(rep_sactimes(2,sac)) - rep_x(rep_sactimes(1,sac));
            sacy = rep_y(rep_sactimes(2,sac)) - rep_y(rep_sactimes(1,sac));
            sacamps(2*img,sac) = sqrt(sacx^2+sacy^2);
        end
    end
    
    %---Store data across all sets---%
    for type = 1:4
        if sum(type == trial_types) > 1
            %post_pupil{type} = [post_pupil{type}; nanmean(pupil_data(trial_types == type,:))];
            post_sacamp{type} = [post_sacamp{type}; nanmean(sacamps(trial_types == type,:))];
            post_fixdur{type} = [post_fixdur{type}; nanmean(fixdurs(trial_types == type,:))];
        end
    end
end
post_time_sets = post_time_sets/post_num_sets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---PLOTS for ROI Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clrs = 'brmgg';

figure
hold on
for r = [1:3 5]
    plot(100*filtfilt(1/200*ones(1,200),1,pre_time_sets(r,:)),...
        [clrs(r) '-']);
end

for r = [1:3 5]
    plot(100*filtfilt(1/200*ones(1,200),1,post_time_sets(r,:)),...
        [clrs(r) '--']);
end
ylim([0 60])
legend('Pre-Novel','Pre-Repeat','Pre-Replaced','Pre-Moved','Post-Novel','Post-Repeat','Post-Replaced','Post-Moved');
title([pre_files{1}(1:2) ' SCMRM Pre vs. Post: % of Time in ROI over Time'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_Percent_Time'])


figure
subplot(2,3,1)
hold on
bar([mean(pre_sets_10numFix);mean(post_sets_10numFix)]','grouped');
errorb([mean(pre_sets_10numFix);mean(post_sets_10numFix)]',...
    [std(pre_sets_10numFix)./sqrt(pre_num_sets); std(post_sets_10numFix)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_10numFix(:,i),post_sets_10numFix(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_10numFix(:,i))+std(pre_sets_10numFix(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of 1st 10 Fixations in ROI')

subplot(2,3,2)
hold on
bar([mean(pre_sets_Time3);mean(post_sets_Time3)]','grouped');
errorb([mean(pre_sets_Time3);mean(post_sets_Time3)]',...
    [std(pre_sets_Time3)./sqrt(pre_num_sets); std(post_sets_Time3)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time3(:,i),post_sets_Time3(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time3(:,i))+std(pre_sets_Time3(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of Time in 1st 3 seconds in ROI')

subplot(2,3,3)
hold on
bar([mean(pre_sets_Time15);mean(post_sets_Time15)]','grouped');
errorb([mean(pre_sets_Time15);mean(post_sets_Time15)]',...
    [std(pre_sets_Time15)./sqrt(pre_num_sets); std(post_sets_Time15)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time15(:,i),post_sets_Time15(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time15(:,i))+std(pre_sets_Time15(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of Time in 0.5-1.5 sec window in ROI')

subplot(2,3,4)
hold on
bar([mean(pre_sets_numFix);mean(post_sets_numFix)]','grouped');
errorb([mean(pre_sets_numFix);mean(post_sets_numFix)]',...
    [std(pre_sets_numFix)./sqrt(pre_num_sets); std(post_sets_numFix)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_numFix(:,i),post_sets_numFix(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_numFix(:,i))+std(pre_sets_numFix(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Fixations in ROI')

subplot(2,3,5)
hold on
bar([mean(pre_sets_Time);mean(post_sets_Time)]','grouped');
errorb([mean(pre_sets_Time);mean(post_sets_Time)]',...
    [std(pre_sets_Time)./sqrt(pre_num_sets); std(post_sets_Time)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time(:,i),post_sets_Time(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time(:,i))+std(pre_sets_Time(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Percentage')
title('% of All Time in ROI')

subtitle([pre_files{1}(1:2) ': Absolute Pre vs. Post SCMRM'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_Absolute_PrePost'])


%normalize by repeat condition
pre_sets_10numFix = pre_sets_10numFix./mean(pre_sets_10numFix(:,1));
post_sets_10numFix = post_sets_10numFix./mean(post_sets_10numFix(:,1));
pre_sets_numFix = pre_sets_numFix./mean(pre_sets_numFix(:,1));
post_sets_numFix = post_sets_numFix./mean(post_sets_numFix(:,1));
pre_sets_Time3 = pre_sets_Time3./mean(pre_sets_Time3(:,1));
post_sets_Time3 = post_sets_Time3./mean(post_sets_Time3(:,1));
pre_sets_Time15 = pre_sets_Time15./mean(pre_sets_Time15(:,1));
post_sets_Time15 = post_sets_Time15./mean(post_sets_Time15(:,1));
pre_sets_Time = pre_sets_Time./mean(pre_sets_Time(:,1));
post_sets_Time = post_sets_Time./mean(post_sets_Time(:,1));

figure
subplot(2,3,1)
hold on
bar([mean(pre_sets_10numFix);mean(post_sets_10numFix)]','grouped');
errorb([mean(pre_sets_10numFix);mean(post_sets_10numFix)]',...
    [std(pre_sets_10numFix)./sqrt(pre_num_sets); std(post_sets_10numFix)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_10numFix(:,i),post_sets_10numFix(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_10numFix(:,i))+std(pre_sets_10numFix(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Relative Percentage to Novel Condition')
title('% of 1st 10 Fixations in ROI')

subplot(2,3,2)
hold on
bar([mean(pre_sets_Time3);mean(post_sets_Time3)]','grouped');
errorb([mean(pre_sets_Time3);mean(post_sets_Time3)]',...
    [std(pre_sets_Time3)./sqrt(pre_num_sets); std(post_sets_Time3)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time3(:,i),post_sets_Time3(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time3(:,i))+std(pre_sets_Time3(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Relative Percentage to Novel Condition')
title('% of Time in 1st 3 seconds in ROI')

subplot(2,3,3)
hold on
bar([mean(pre_sets_Time15);mean(post_sets_Time15)]','grouped');
errorb([mean(pre_sets_Time15);mean(post_sets_Time15)]',...
    [std(pre_sets_Time15)./sqrt(pre_num_sets); std(post_sets_Time15)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time15(:,i),post_sets_Time15(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time15(:,i))+std(pre_sets_Time15(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Relative Percentage to Novel Condition')
title('% of Time in 0.5-1.5 sec window in ROI')

subplot(2,3,4)
hold on
bar([mean(pre_sets_numFix);mean(post_sets_numFix)]','grouped');
errorb([mean(pre_sets_numFix);mean(post_sets_numFix)]',...
    [std(pre_sets_numFix)./sqrt(pre_num_sets); std(post_sets_numFix)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_numFix(:,i),post_sets_numFix(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_numFix(:,i))+std(pre_sets_numFix(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Relative Percentage to Novel Condition')
title('% of All Fixations in ROI')

subplot(2,3,5)
hold on
bar([mean(pre_sets_Time);mean(post_sets_Time)]','grouped');
errorb([mean(pre_sets_Time);mean(post_sets_Time)]',...
    [std(pre_sets_Time)./sqrt(pre_num_sets); std(post_sets_Time)./sqrt(post_num_sets)]')
for i = 1:4
    [~,p] = ttest2(pre_sets_Time(:,i),post_sets_Time(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time(:,i))+std(pre_sets_Time(:,i)),'k*')
    end
end
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Relative Percentage to Novel Condition')
title('% of All Time in ROI')

subtitle([pre_files{1}(1:2) ': Relative Pre vs. Post SCMRM'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_Relative_PrePost'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---PLOTS for Eye Movements---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Absolute Fixation Duration---%
figure
subplot(2,2,1)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{2}(:,1:20)),nanstd(pre_fixdur{2}(:,1:20))./sqrt(size(pre_fixdur{2},1)),'r')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{2}(:,1:20)),nanstd(post_fixdur{2}(:,1:20))./sqrt(size(post_fixdur{2},1)),'--r')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Repeat')

subplot(2,2,2)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{3}(:,1:20)),nanstd(pre_fixdur{3}(:,1:20))./sqrt(size(pre_fixdur{3},1)),'m')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{3}(:,1:20)),nanstd(post_fixdur{3}(:,1:20))./sqrt(size(post_fixdur{3},1)),'--m')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Replaced')

subplot(2,2,3)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{4}(:,1:20)),nanstd(pre_fixdur{4}(:,1:20))./sqrt(size(pre_fixdur{4},1)),'g')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{4}(:,1:20)),nanstd(post_fixdur{4}(:,1:20))./sqrt(size(post_fixdur{4},1)),'--g')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Moved')

pre_2nd = [pre_fixdur{2} pre_fixdur{3} pre_fixdur{4}];
post_2nd = [post_fixdur{2} post_fixdur{3} post_fixdur{4}];

subplot(2,2,4)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(pre_2nd(:,1:20)),nanstd(pre_2nd(:,1:20))./sqrt(size(pre_fixdur{1},1)),'k')
errorbar(nanmean(post_2nd(:,1:20)),nanstd(post_2nd(:,1:20))./sqrt(size(post_fixdur{1},1)),'--k')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs 2nd Presentation')

subtitle([pre_files{1}(1:2) 'SCMRM Pre vs Post: Fixation Duration'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_FixationDuration'])


%---Relative Fixation Duration---%
%normalize to First Fixation on Novel images
for c = 4:-1:1
    pre_fixdur{c}  = pre_fixdur{c}./nanmean(pre_fixdur{1}(:,1));
    post_fixdur{c} = post_fixdur{c}./nanmean(post_fixdur{1}(:,1));
end

figure
subplot(2,2,1)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{2}(:,1:20)),nanstd(pre_fixdur{2}(:,1:20))./sqrt(size(pre_fixdur{2},1)),'r')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{2}(:,1:20)),nanstd(post_fixdur{2}(:,1:20))./sqrt(size(post_fixdur{2},1)),'--r')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Repeat')

subplot(2,2,2)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{3}(:,1:20)),nanstd(pre_fixdur{3}(:,1:20))./sqrt(size(pre_fixdur{3},1)),'m')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{3}(:,1:20)),nanstd(post_fixdur{3}(:,1:20))./sqrt(size(post_fixdur{3},1)),'--m')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Replaced')

subplot(2,2,3)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(pre_fixdur{4}(:,1:20)),nanstd(pre_fixdur{4}(:,1:20))./sqrt(size(pre_fixdur{4},1)),'g')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(post_fixdur{4}(:,1:20)),nanstd(post_fixdur{4}(:,1:20))./sqrt(size(post_fixdur{4},1)),'--g')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs Moved')

pre_2nd = [pre_fixdur{2} pre_fixdur{3} pre_fixdur{4}];
post_2nd = [post_fixdur{2} post_fixdur{3} post_fixdur{4}];

subplot(2,2,4)
hold on
errorbar(nanmean(pre_fixdur{1}(:,1:20)),nanstd(pre_fixdur{1}(:,1:20))./sqrt(size(pre_fixdur{1},1)),'b')
errorbar(nanmean(post_fixdur{1}(:,1:20)),nanstd(post_fixdur{1}(:,1:20))./sqrt(size(post_fixdur{1},1)),'--b')
errorbar(nanmean(pre_2nd(:,1:20)),nanstd(pre_2nd(:,1:20))./sqrt(size(pre_fixdur{1},1)),'k')
errorbar(nanmean(post_2nd(:,1:20)),nanstd(post_2nd(:,1:20))./sqrt(size(post_fixdur{1},1)),'--k')
hold off
xlim([0 20.5])
xlabel('Ordinal Fixation Number')
ylabel('Fixation Duration (ms)')
title('Novel vs 2nd Presentation')

subtitle([pre_files{1}(1:2) 'SCMRM Pre vs Post: Relative Fixation Duration'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_Relative_FixationDuration'])


%---Absolute Saccade Amplitudes---%
figure
subplot(2,2,1)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20))/24,nanstd(pre_sacamp{1}(:,1:20))/24./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{2}(:,1:20))/24,nanstd(pre_sacamp{2}(:,1:20))/24./sqrt(size(pre_sacamp{2},1)),'r')
errorbar(nanmean(post_sacamp{1}(:,1:20))/24,nanstd(post_sacamp{1}(:,1:20))/24./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{2}(:,1:20))/24,nanstd(post_sacamp{2}(:,1:20))/24./sqrt(size(post_sacamp{2},1)),'--r')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Ampitude (dva)')
title('Novel vs Repeat')

subplot(2,2,2)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20))/24,nanstd(pre_sacamp{1}(:,1:20))/24./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{3}(:,1:20))/24,nanstd(pre_sacamp{3}(:,1:20))/24./sqrt(size(pre_sacamp{3},1)),'m')
errorbar(nanmean(post_sacamp{1}(:,1:20))/24,nanstd(post_sacamp{1}(:,1:20))/24./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{3}(:,1:20))/24,nanstd(post_sacamp{3}(:,1:20))/24./sqrt(size(post_sacamp{3},1)),'--m')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Ampitude (dva)')
title('Novel vs Replaced')

subplot(2,2,3)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20))/24,nanstd(pre_sacamp{1}(:,1:20))/24./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{4}(:,1:20))/24,nanstd(pre_sacamp{4}(:,1:20))/24./sqrt(size(pre_sacamp{4},1)),'g')
errorbar(nanmean(post_sacamp{1}(:,1:20))/24,nanstd(post_sacamp{1}(:,1:20))/24./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{4}(:,1:20))/24,nanstd(post_sacamp{4}(:,1:20))/24./sqrt(size(post_sacamp{4},1)),'--g')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Ampitude (dva)')
title('Novel vs Moved')

pre_2nd = [pre_sacamp{2} pre_sacamp{3} pre_sacamp{4}];
post_2nd = [post_sacamp{2} post_sacamp{3} post_sacamp{4}];

subplot(2,2,4)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20))/24,nanstd(pre_sacamp{1}(:,1:20))/24./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(post_sacamp{1}(:,1:20))/24,nanstd(post_sacamp{1}(:,1:20))/24./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(pre_2nd(:,1:20))/24,nanstd(pre_2nd(:,1:20))/24./sqrt(size(pre_sacamp{1},1)),'k')
errorbar(nanmean(post_2nd(:,1:20))/24,nanstd(post_2nd(:,1:20))/24./sqrt(size(post_sacamp{1},1)),'--k')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Ampitude (dva)')
title('Novel vs 2nd Presentation')

subtitle([pre_files{1}(1:2) ' SCMRM Pre vs Post: Saccade Amplitude'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_SaccadeAmplitude'])

%---Relative Fixation Duration---%
%normalize to First Saccade on Novel images
for c = 4:-1:1
    pre_sacamp{c}  = pre_sacamp{c}./nanmean(pre_sacamp{1}(:,1));
    post_sacamp{c} = post_sacamp{c}./nanmean(post_sacamp{1}(:,1));
end

figure
subplot(2,2,1)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20)),nanstd(pre_sacamp{1}(:,1:20))./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{2}(:,1:20)),nanstd(pre_sacamp{2}(:,1:20))./sqrt(size(pre_sacamp{2},1)),'r')
errorbar(nanmean(post_sacamp{1}(:,1:20)),nanstd(post_sacamp{1}(:,1:20))./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{2}(:,1:20)),nanstd(post_sacamp{2}(:,1:20))./sqrt(size(post_sacamp{2},1)),'--r')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Amplitude (dva)')
title('Novel vs Repeat')

subplot(2,2,2)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20)),nanstd(pre_sacamp{1}(:,1:20))./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{3}(:,1:20)),nanstd(pre_sacamp{3}(:,1:20))./sqrt(size(pre_sacamp{3},1)),'m')
errorbar(nanmean(post_sacamp{1}(:,1:20)),nanstd(post_sacamp{1}(:,1:20))./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{3}(:,1:20)),nanstd(post_sacamp{3}(:,1:20))./sqrt(size(post_sacamp{3},1)),'--m')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Amplitude (dva)')
title('Novel vs Replaced')

subplot(2,2,3)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20)),nanstd(pre_sacamp{1}(:,1:20))./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(pre_sacamp{4}(:,1:20)),nanstd(pre_sacamp{4}(:,1:20))./sqrt(size(pre_sacamp{4},1)),'g')
errorbar(nanmean(post_sacamp{1}(:,1:20)),nanstd(post_sacamp{1}(:,1:20))./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(post_sacamp{4}(:,1:20)),nanstd(post_sacamp{4}(:,1:20))./sqrt(size(post_sacamp{4},1)),'--g')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Amplitude (dva)')
title('Novel vs Moved')

pre_2nd = [pre_sacamp{2} pre_sacamp{3} pre_sacamp{4}];
post_2nd = [post_sacamp{2} post_sacamp{3} post_sacamp{4}];

subplot(2,2,4)
hold on
errorbar(nanmean(pre_sacamp{1}(:,1:20)),nanstd(pre_sacamp{1}(:,1:20))./sqrt(size(pre_sacamp{1},1)),'b')
errorbar(nanmean(post_sacamp{1}(:,1:20)),nanstd(post_sacamp{1}(:,1:20))./sqrt(size(post_sacamp{1},1)),'--b')
errorbar(nanmean(pre_2nd(:,1:20)),nanstd(pre_2nd(:,1:20))./sqrt(size(pre_sacamp{1},1)),'k')
errorbar(nanmean(post_2nd(:,1:20)),nanstd(post_2nd(:,1:20))./sqrt(size(post_sacamp{1},1)),'--k')
hold off
xlim([0 20.5])
xlabel('Ordinal Saccade Number')
ylabel('Saccade Amplitude (dva)')
title('Novel vs 2nd Presentation')

subtitle([pre_files{1}(1:2) 'SCMRM Pre vs Post: Relative Saccade Amplitude'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_Relative_SaccadeAmplitude'])


%---Plot Pupil Diameter Over Time---%
t = 1:5:7500;
figure
subplot(1,2,1)
hold on
for c = 1:4
    dofill(t,pre_pupil{c}/1000,clrs(c),1,60)
end
yl = ylim;
plot([500 500],[yl(1) yl(2)],'--k')
hold off
ylim(yl)
xlim([0 7500])
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
title('Pre-Lesion')

subplot(1,2,2)
hold on
for c = 1:4
    dofill(t,post_pupil{c}/1000,clrs(c),1,60)
end
yl = ylim;
plot([500 500],[yl(1) yl(2)],'--k')
hold off
xlim([0 7500])
ylim(yl)
xlabel('Time from Image Onset (ms)')
ylabel('Normalize Pupil Data (a.u.)')
set(gca,'Xtick',0:500:7500)
set(gca,'XtickLabel',num2cell([0:500:7500]-500))
title('Post-Lesion')
legend('Novel','Repeat','Replaced','Moved')