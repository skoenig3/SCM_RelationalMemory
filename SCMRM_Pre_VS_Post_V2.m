%compare pre-post lesion
%orginal written by Seth Konig 6/7/16 updated to V2 on10/17/2017
clar
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

 
%---Tobii---%
% pre_files = {'TO150901.2','TO150903.2','TO150904.3',...
%     'TO150908.2','TO150909.2','TO150910.2','TO150911.2',...
%     'TO150914.2','TO150915.2','TO150916.2','TO150918.2',...
%     'TO151001.2','TO151002.2','TO151006.2','TO151007.2'};
% 
% post_files = {'TO170620.2','TO170621.2','TO170622.2','TO170623.2',...
%     'TO170626.2','TO170627.2','TO170628.2','TO170629.2',...
%     'TO170630.2','TO170705.2','TO170706.2','TO170707.2',...
%     'TO170710.2','TO170711.2','TO170712.2'};

%---Manfred---%
pre_files = {'MF170222.2','MF170223.2','MF170224.2','MF170227.2',...
    'MF170301.2','MF170303.2','MF170306.2','MF170307.2',...
    'MF170308.2','MF170309.2','MF170313.2','MF170314.2',... 
    'MF170315.2','MF170316.2','MF170317.2'};
post_files = {'MF170222.2','MF170223.2','MF170224.2','MF170227.2',...
    'MF170301.2','MF170303.2','MF170306.2','MF170307.2',...
    'MF170308.2','MF170309.2','MF170313.2','MF170314.2',...
    'MF170315.2','MF170316.2','MF170317.2'};

%---All Combined---%
% pre_files = {'TO150901.2','TO150903.2','TO150904.3',...
%     'TO150908.2','TO150909.2','TO150910.2','TO150911.2',...
%     'TO150914.2','TO150915.2','TO150916.2','TO150918.2',...
%     'TO151001.2','TO151002.2','TO151006.2','TO151007.2',...
%     'RR151027.2','RR151026.2','RR151023.2','RR151022.2',...
%     'RR151021.2','RR151020.2','RR151019.2','RR151016.2',...
%     'RR151015.3','RR151014.2','RR151013.2','RR151009.2',...
%     'RR151008.2','RR151007.2','RR151005.2','RR151002.2',...
%     'RR151001.2',...
%     'PW150901.2','PW150902.2','PW150903.2','PW150904.2',...
%     'PW150908.2','PW150910.2','PW150911.2',...
%     'PW150914.2','PW150916.2','PW150917.2','PW150918.2',...
%     'PW150921.2','PW150922.2','PW150923.2','PW150924.2'};
% 
% post_files = {'TO170620.2','TO170621.2','TO170622.2','TO170623.2',...
%     'TO170626.2','TO170627.2','TO170628.2','TO170629.2',...
%     'TO170630.2','TO170705.2','TO170706.2','TO170707.2',...
%     'TO170710.2','TO170711.2','TO170712.2',...
%     'RR160425.2','RR160426.2','RR160701.2','RR160706.2','RR160707.2',...
%     'RR160708.2','RR160714.2','RR160715.2','RR160718.2','RR160721.2',...
%     'RR160725.2','RR160726.2','RR160727.2','RR160728.2','RR160729.2',...
%     'RR160801.2','RR160803.2','RR160804.2','RR160805.2',...
%     'PW160401.2','PW160405.2','PW160406.2','PW160407.2','PW160408.2',...
%     'PW160411.2','PW160412.2','PW160413.2','PW160414.2','PW160415.2',...
%     'PW160418.2','PW160419.2','PW160420.2','PW160421.2','PW160422.2'};


%---Get Pre-Data---%
pre_num_sets = length(pre_files);
pre_sets_numFix = NaN(length(pre_files),4);
pre_sets_10numFix = NaN(length(pre_files),4);
pre_sets_Time = NaN(length(pre_files),4);
pre_sets_Time3 = NaN(length(pre_files),4);
pre_sets_Time15 = NaN(length(pre_files),4);
pre_time_sets = zeros(5,7000);

%Eye data
pre_pupil = cell(1,2);
pre_sacamp = cell(1,2);
pre_sacamp{1} = NaN(length(pre_files),50);
pre_sacamp{2} = NaN(length(pre_files),50);
pre_fixdur = cell(1,2);
pre_fixdur{1} = NaN(length(pre_files),50);
pre_fixdur{2} = NaN(length(pre_files),50);

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
        
        if nov_img_off-nov_img_on > 10500
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
        
        if rep_img_off-rep_img_on > 10050
            continue
        end
        
        rep_x = fixationstats{2*img}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{2*img}.XY(2,rep_img_on:rep_img_off);
        try %may happen do to EOG overflow?
            rep_pupil = pupildata{2*img}(round(rep_img_on/5):round(rep_img_off/5));
        catch
            rep_pupil = NaN(1,round(rep_img_off/5)-round(rep_img_on/5));
        end
        
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
        try %may happen do to EOG overflow?
            rep_pupil = pupildata{2*img}(round(rep_img_on/5)-100:round(rep_img_off/5));
        catch
            rep_pupil = NaN(1,round(rep_img_off/5)-100-round(rep_img_on/5));
        end
       
        
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
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil+2000;
        rep_pupil = rep_pupil+2000;
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)));
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)));
        
        %only look at the 1st 7 seconds the images is up
        nov_pupil = nov_pupil(1:1500);
        try %may happen do to EOG overflow?
            rep_pupil = rep_pupil(1:1500);
        catch
            rep_pupil = [rep_pupil; NaN(1500-length(rep_pupil),1)];
        end
        
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
    for type = 1:2
        if type == 1
            median_fix_count = floor(median(sum(~isnan(fixdurs(trial_types == type,:))')));
            median_sac_count = floor(median(sum(~isnan(sacamps(trial_types == type,:))')));
            pre_pupil{type} = [pre_pupil{type}; nanmean(pupil_data(trial_types == type,:))];
            pre_sacamp{type}(sets,1:median_sac_count) = nanmean(sacamps(trial_types == type,1:median_sac_count));
            pre_fixdur{type}(sets,1:median_fix_count) = nanmean(fixdurs(trial_types == type,1:median_fix_count));
        else %all repeats
            median_fix_count = floor(median(sum(~isnan(fixdurs(trial_types >= type,:))')));
            median_sac_count = floor(median(sum(~isnan(sacamps(trial_types >= type,:))')));
            pre_pupil{type} = [pre_pupil{type}; nanmean(pupil_data(trial_types >= type,:))];
            pre_sacamp{type}(sets,1:median_sac_count) = nanmean(sacamps(trial_types >= type,1:median_sac_count));
            pre_fixdur{type}(sets,1:median_fix_count) = nanmean(fixdurs(trial_types >= type,1:median_fix_count));
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
post_pupil = cell(1,2);
post_sacamp = cell(1,2);
post_sacamp{1} = NaN(length(post_files),50);
post_sacamp{2} = NaN(length(post_files),50);
post_fixdur = cell(1,2);
post_fixdur{1} = NaN(length(post_files),50);
post_fixdur{2} = NaN(length(post_files),50);
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
        
        if nov_img_off-nov_img_on > 10500
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
        
        if rep_img_off-rep_img_on > 10050
            continue
        end
        
        rep_x = fixationstats{2*img}.XY(1,rep_img_on:rep_img_off);
        rep_y = fixationstats{2*img}.XY(2,rep_img_on:rep_img_off);
        try %may happen do to EOG overflow?
            rep_pupil = pupildata{2*img}(round(rep_img_on/5):round(rep_img_off/5));
        catch
            rep_pupil = NaN(1,round(rep_img_off/5)-round(rep_img_on/5));
        end
        
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
        try %may happen do to EOG overflow?
            rep_pupil = pupildata{2*img}(round(rep_img_on/5)-100:round(rep_img_off/5));
        catch
            rep_pupil = NaN(1,round(rep_img_off/5)-100-round(rep_img_on/5));
        end
       
        
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
        
        %normalize to 100 ms of prestimulus levels
        nov_pupil = nov_pupil+2000;
        rep_pupil = rep_pupil+2000;
        nov_pupil = nov_pupil./abs(mean(nov_pupil(81:100)));
        rep_pupil = rep_pupil./abs(mean(rep_pupil(81:100)));
        
        %only look at the 1st 7 seconds the images is up
        nov_pupil = nov_pupil(1:1500);
        nov_pupil = nov_pupil(1:1500);
        try %may happen do to EOG overflow?
            rep_pupil = rep_pupil(1:1500);
        catch
            rep_pupil = [rep_pupil; NaN(1500-length(rep_pupil),1)];
        end
        
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
    for type = 1:2
        if type == 1
            median_fix_count = floor(median(sum(~isnan(fixdurs(trial_types == type,:))')));
            median_sac_count = floor(median(sum(~isnan(sacamps(trial_types == type,:))')));
            post_pupil{type} = [post_pupil{type}; nanmean(pupil_data(trial_types == type,:))];
            post_sacamp{type}(sets,1:median_sac_count) = nanmean(sacamps(trial_types == type,1:median_sac_count));
            post_fixdur{type}(sets,1:median_fix_count) = nanmean(fixdurs(trial_types == type,1:median_fix_count));
        else %all repeats
            median_fix_count = floor(median(sum(~isnan(fixdurs(trial_types >= type,:))')));
            median_sac_count = floor(median(sum(~isnan(sacamps(trial_types >= type,:))')));
            post_pupil{type} = [post_pupil{type}; nanmean(pupil_data(trial_types >= type,:))];
            post_sacamp{type}(sets,1:median_sac_count) = nanmean(sacamps(trial_types >= type,1:median_sac_count));
            post_fixdur{type}(sets,1:median_fix_count) = nanmean(fixdurs(trial_types >= type,1:median_fix_count));
        end
    end
end
post_time_sets = post_time_sets/post_num_sets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---PLOTS for ROI Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clrs = 'brmgg';
%%
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
ylabel('% of Time Looking in ROI')
legend('Pre-Novel','Pre-Repeat','Pre-Replaced','Pre-Moved','Post-Novel','Post-Repeat','Post-Replaced','Post-Moved');

title([pre_files{1}(1:2) ' SCMRM Pre vs. Post: Time Course'])
%save_and_close_fig(figure_dir,[pre_files{1}(1:2) '_SCMRM_PrePost_Percent_Time'])

%%

npre = size(pre_sets_Time15,1);
npost = size(post_sets_Time15,1);
ids = [[ones(npre,1); 2*ones(npre,1); 3*ones(npre,1); 4*ones(npre,1); ...
    ones(npost,1); 2*ones(npost,1); 3*ones(npost,1); 4*ones(npost,1)]...
    [ones(npre,1); ones(npre,1); ones(npre,1); ones(npre,1);...
    2*ones(npost,1); 2*ones(npost,1); 2*ones(npost,1); 2*ones(npost,1)]];

vals = [pre_sets_Time15(:,1); pre_sets_Time15(:,2); pre_sets_Time15(:,3); pre_sets_Time15(:,4);...
    post_sets_Time15(:,1); post_sets_Time15(:,2); post_sets_Time15(:,3); post_sets_Time15(:,4)];
[P_ANOVA_T15] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [pre_sets_Time(:,1); pre_sets_Time(:,2); pre_sets_Time(:,3); pre_sets_Time(:,4);...
    post_sets_Time(:,1); post_sets_Time(:,2); post_sets_Time(:,3); post_sets_Time(:,4)];
[P_ANOVA_T_all] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [pre_sets_numFix(:,1); pre_sets_numFix(:,2); pre_sets_numFix(:,3); pre_sets_numFix(:,4);...
    post_sets_numFix(:,1); post_sets_numFix(:,2); post_sets_numFix(:,3); post_sets_numFix(:,4)];
[P_ANOVA_F_all] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures
%%
figure
subplot(1,3,1)
errorb([mean(pre_sets_Time15);mean(post_sets_Time15)]',...
    [std(pre_sets_Time15)./sqrt(pre_num_sets); std(post_sets_Time15)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(pre_sets_Time15(:,i),post_sets_Time15(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time15(:,i))+std(pre_sets_Time15(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of Time from 0.5-1.5 sec')
title(sprintf(['Time 0.5-1.5 seeconds \n p_{lesion} = ' num2str(P_ANOVA_T15(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_T15(1),2) ', p_{inter} = ' num2str(P_ANOVA_T15(3),2)]))

subplot(1,3,2)
errorb([mean(pre_sets_numFix);mean(post_sets_numFix)]',...
    [std(pre_sets_numFix)./sqrt(pre_num_sets); std(post_sets_numFix)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(pre_sets_numFix(:,i),post_sets_numFix(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_numFix(:,i))+std(pre_sets_numFix(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4) 
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of All Fixations ')
title(sprintf(['Percent of All Fixations \n p_{lesion} = ' num2str(P_ANOVA_F_all(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_F_all(1),2) ', p_{inter} = ' num2str(P_ANOVA_F_all(3),2)]))

subplot(1,3,3)
errorb([mean(pre_sets_Time);mean(post_sets_Time)]',...
    [std(pre_sets_Time)./sqrt(pre_num_sets); std(post_sets_Time)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(pre_sets_Time(:,i),post_sets_Time(:,i));
    if p < 0.05
        plot(i, mean(pre_sets_Time(:,i))+std(pre_sets_Time(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of All Time')
title(sprintf(['Percent of All Time \n p_{lesion} = ' num2str(P_ANOVA_T_all(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_T_all(1),2) ', p_{inter} = ' num2str(P_ANOVA_T_all(3),2)]))

%subtitle([pre_files{1}(1:2) ': Raw Pre vs. Post SCMRM'])

%%
%normalize by to Novel condition
norm_pre_sets_Time15 = pre_sets_Time15./mean(pre_sets_Time15(:,1));
norm_pre_sets_numFix = pre_sets_numFix./mean(pre_sets_numFix(:,1));
norm_pre_sets_Time = pre_sets_Time./mean(pre_sets_Time(:,1));
norm_post_sets_Time15 = post_sets_Time15./mean(post_sets_Time15(:,1));
norm_post_sets_numFix = post_sets_numFix./mean(post_sets_numFix(:,1));
norm_post_sets_Time = post_sets_Time./mean(post_sets_Time(:,1));


vals = [norm_pre_sets_Time15(:,1); norm_pre_sets_Time15(:,2); norm_pre_sets_Time15(:,3); norm_pre_sets_Time15(:,4);...
    norm_post_sets_Time15(:,1); norm_post_sets_Time15(:,2); norm_post_sets_Time15(:,3); norm_post_sets_Time15(:,4)];
[P_ANOVA_T15] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [norm_pre_sets_Time(:,1); norm_pre_sets_Time(:,2); norm_pre_sets_Time(:,3); norm_pre_sets_Time(:,4);...
    norm_post_sets_Time(:,1); norm_post_sets_Time(:,2); norm_post_sets_Time(:,3); norm_post_sets_Time(:,4)];
[P_ANOVA_T_all] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

vals = [norm_pre_sets_numFix(:,1); norm_pre_sets_numFix(:,2); norm_pre_sets_numFix(:,3); norm_pre_sets_numFix(:,4);...
    norm_post_sets_numFix(:,1); norm_post_sets_numFix(:,2); norm_post_sets_numFix(:,3); norm_post_sets_numFix(:,4)];
[P_ANOVA_F_all] = anovan(vals,ids,'model','interaction','varnames',{'Type','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

figure
subplot(1,3,1)
errorb([mean(norm_pre_sets_Time15);mean(norm_post_sets_Time15)]',...
    [std(norm_pre_sets_Time15)./sqrt(pre_num_sets); std(norm_post_sets_Time15)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(norm_pre_sets_Time15(:,i),norm_post_sets_Time15(:,i));
    if p < 0.05
        plot(i, mean(norm_pre_sets_Time15(:,i))+std(norm_pre_sets_Time15(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of Time from 0.5-1.5 sec')
title(sprintf(['Time 0.5-1.5 seeconds \n p_{lesion} = ' num2str(P_ANOVA_T15(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_T15(1),2) ', p_{inter} = ' num2str(P_ANOVA_T15(3),2)]))

subplot(1,3,2)
errorb([mean(norm_pre_sets_numFix);mean(norm_post_sets_numFix)]',...
    [std(norm_pre_sets_numFix)./sqrt(pre_num_sets); std(norm_post_sets_numFix)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(norm_pre_sets_numFix(:,i),norm_post_sets_numFix(:,i));
    if p < 0.05
        plot(i, mean(norm_pre_sets_numFix(:,i))+std(norm_pre_sets_numFix(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of All Fixations ')
title(sprintf(['Percent of All Fixations \n p_{lesion} = ' num2str(P_ANOVA_F_all(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_F_all(1),2) ', p_{inter} = ' num2str(P_ANOVA_F_all(3),2)]))

subplot(1,3,3)
errorb([mean(norm_pre_sets_Time);mean(norm_post_sets_Time)]',...
    [std(norm_pre_sets_Time)./sqrt(pre_num_sets); std(norm_post_sets_Time)./sqrt(post_num_sets)]')
hold on
for i = 1:4
    [~,p] = ttest2(norm_pre_sets_Time(:,i),norm_post_sets_Time(:,i));
    if p < 0.05
        plot(i, mean(norm_pre_sets_Time(:,i))+std(norm_pre_sets_Time(:,i)),'k*')
    end
end
hold off
box off
axis square
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('% of All Time')
title(sprintf(['Percent of All Time \n p_{lesion} = ' num2str(P_ANOVA_T_all(2),2) ...
    ' , p_{image type} = ' num2str(P_ANOVA_T_all(1),2) ', p_{inter} = ' num2str(P_ANOVA_T_all(3),2)]))

subtitle([pre_files{1}(1:2) ': Normalized to Novel Condition Pre vs. Post SCMRM'])

%%
median_pre_nov_fix_count = median(sum(~isnan(pre_fixdur{1}')));
median_pre_rep_fix_count = median(sum(~isnan(pre_fixdur{2}')));
median_post_nov_fix_count = median(sum(~isnan(post_fixdur{1}')));
median_post_rep_fix_count = median(sum(~isnan(post_fixdur{2}')));


npre = size(pre_fixdur{1},1);
npost = size(post_fixdur{1},1);

pre_nov = mean(pre_fixdur{1}(:,3:15),2);
pre_rep = mean(pre_fixdur{2}(:,3:15),2);
post_nov = mean(post_fixdur{1}(:,3:15),2);
post_rep = mean(post_fixdur{2}(:,3:15),2);

ids = [[ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)] ...
      [ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)]];

vals = [pre_nov; pre_rep; post_nov; post_rep];
[P_ANOVA_fixdurs] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

figure

subplot(2,2,1)
hold all
errorbar(nanmean(pre_fixdur{1}(:,1:median_pre_nov_fix_count)),nanstd(pre_fixdur{1}(:,1:median_pre_nov_fix_count))./sqrt(npre))
errorbar(nanmean(pre_fixdur{2}(:,1:median_pre_rep_fix_count)),nanstd(pre_fixdur{2}(:,1:median_pre_rep_fix_count))./sqrt(npre))
errorbar(nanmean(post_fixdur{1}(:,1:median_post_nov_fix_count)),nanstd(post_fixdur{1}(:,1:median_post_nov_fix_count))./sqrt(npost))
errorbar(nanmean(post_fixdur{2}(:,1:median_post_rep_fix_count)),nanstd(post_fixdur{2}(:,1:median_post_rep_fix_count))./sqrt(npost))
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_fix_count+1])
title(' Fixation Durations')

[~,p_nov] = ttest2(pre_nov,post_nov);
[~,p_rep] = ttest2(pre_rep,post_rep);

subplot(2,2,3)
errorb([[mean(pre_nov);mean(post_nov)] [mean(pre_rep);mean(post_rep)]]',...
    [[std(pre_nov)./sqrt(npre) ;std(post_nov)./sqrt(npost)] ...
    [std(pre_rep)./sqrt(npre) ;std(post_rep)./sqrt(npost)]]')
hold on
if p_nov < 0.05
    plot(1,mean(pre_nov)+2*std(pre_nov),'k*')
end
if p_rep < 0.05
     plot(2,mean(pre_rep)+2*std(pre_rep),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Fixation Duration (ms)')
title(sprintf(['ANOVA Ordinal Fixations 3-15: \n p_{lesion} = ' num2str(P_ANOVA_fixdurs(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_fixdurs(1),2) ', p_{inter} = ' num2str(P_ANOVA_fixdurs(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')


norm_pre_nov_fixdur = pre_fixdur{1}(:,1:median_pre_nov_fix_count);
norm_pre_rep_fixdur = pre_fixdur{2}(:,1:median_pre_rep_fix_count);
norm_pre_rep_fixdur = norm_pre_rep_fixdur./nanmean(norm_pre_nov_fixdur(:,1));
norm_pre_nov_fixdur = norm_pre_nov_fixdur./nanmean(norm_pre_nov_fixdur(:,1));
norm_post_rep_fixdur = post_fixdur{2}(:,1:median_post_rep_fix_count);
norm_post_nov_fixdur = post_fixdur{1}(:,1:median_post_nov_fix_count);
norm_post_rep_fixdur = norm_post_rep_fixdur./nanmean(norm_post_nov_fixdur(:,1));
norm_post_nov_fixdur = norm_post_nov_fixdur./nanmean(norm_post_nov_fixdur(:,1));


norm_pre_nov = mean(norm_pre_nov_fixdur(:,3:15),2);
norm_pre_rep = mean(norm_pre_rep_fixdur(:,3:15),2);
norm_post_nov = mean(norm_post_nov_fixdur(:,3:15),2);
norm_post_rep = mean(norm_post_rep_fixdur(:,3:15),2);

vals = [norm_pre_nov; norm_pre_rep; norm_post_nov; norm_post_rep];
[P_ANOVA_fixdurs] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

subplot(2,2,2)
hold all
errorbar(nanmean(norm_pre_nov_fixdur),nanstd(norm_pre_nov_fixdur)./sqrt(npre))
errorbar(nanmean(norm_pre_rep_fixdur),nanstd(norm_pre_rep_fixdur)./sqrt(npre))
errorbar(nanmean(norm_post_nov_fixdur),nanstd(norm_post_nov_fixdur)./sqrt(npost))
errorbar(nanmean(norm_post_rep_fixdur),nanstd(norm_post_rep_fixdur)./sqrt(npost))
hold off
xlabel('Ordinal Fixation #')
ylabel('Relative Fixation Duration ')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_fix_count+1])
title('Normalized Fixation Durations')


[~,norm_p_nov] = ttest2(norm_pre_nov,norm_post_nov);
[~,norm_p_rep] = ttest2(norm_pre_rep,norm_post_rep);

subplot(2,2,4)
errorb([[mean(norm_pre_nov);mean(norm_post_nov)] [mean(norm_pre_rep);mean(norm_post_rep)]]',...
    [[std(norm_pre_nov)./sqrt(npre) ;std(norm_post_nov)./sqrt(npost)] ...
    [std(norm_pre_rep)./sqrt(npre) ;std(norm_post_rep)./sqrt(npost)]]')
hold on
if norm_p_nov < 0.05
    plot(1,mean(norm_pre_nov)+2*std(norm_pre_nov),'k*')
end
if norm_p_rep < 0.05
     plot(2,mean(norm_pre_rep)+2*std(norm_pre_rep),'k*')
end
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Relative Fixation Duration ')
title(sprintf(['ANOVA Ordinal Fixations 3-15: \n p_{lesion} = ' num2str(P_ANOVA_fixdurs(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_fixdurs(1),2) ', p_{inter} = ' num2str(P_ANOVA_fixdurs(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')
%%
median_pre_nov_sac_count = median(sum(~isnan(pre_sacamp{1}')));
median_pre_rep_sac_count = median(sum(~isnan(pre_sacamp{2}')));
median_post_nov_sac_count = median(sum(~isnan(post_sacamp{1}')));
median_post_rep_sac_count = median(sum(~isnan(post_sacamp{2}')));


npre = size(pre_sacamp{1},1);
npost = size(post_sacamp{1},1);

pre_nov = nanmean(pre_sacamp{1}(:,2:5),2)/24;
pre_rep = nanmean(pre_sacamp{2}(:,2:5),2)/24;
post_nov = nanmean(post_sacamp{1}(:,2:5),2)/24;
post_rep = nanmean(post_sacamp{2}(:,2:5),2)/24;

ids = [[ones(npre,1); 2*ones(npre,1); ones(npost,1); 2*ones(npost,1)] ...
      [ones(npre,1); ones(npre,1); 2*ones(npost,1); 2*ones(npost,1)]];

vals = [pre_nov; pre_rep; post_nov; post_rep];
[P_ANOVA_sacamps] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

figure

subplot(2,2,1)
hold all
errorbar(nanmean(pre_sacamp{1}(:,1:median_pre_nov_sac_count))/24,nanstd(pre_sacamp{1}(:,1:median_pre_nov_sac_count))/24./sqrt(npre))
errorbar(nanmean(pre_sacamp{2}(:,1:median_pre_rep_sac_count))/24,nanstd(pre_sacamp{2}(:,1:median_pre_rep_sac_count))/24./sqrt(npre))
errorbar(nanmean(post_sacamp{1}(:,1:median_post_nov_sac_count))/24,nanstd(post_sacamp{1}(:,1:median_post_nov_sac_count))/24./sqrt(npost))
errorbar(nanmean(post_sacamp{2}(:,1:median_post_rep_sac_count))/24,nanstd(post_sacamp{2}(:,1:median_post_rep_sac_count))/24./sqrt(npost))
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_sac_count+1])
title(' Saccade Amplitudes')

[~,p_nov] = ttest2(pre_nov,post_nov);
[~,p_rep] = ttest2(pre_rep,post_rep);

subplot(2,2,3)
errorb([[mean(pre_nov);mean(post_nov)] [mean(pre_rep);mean(post_rep)]]',...
    [[std(pre_nov)./sqrt(npre) ;std(post_nov)./sqrt(npost)] ...
    [std(pre_rep)./sqrt(npre) ;std(post_rep)./sqrt(npost)]]')
hold on
if p_nov < 0.05
    plot(1,mean(pre_nov)+2*std(pre_nov),'k*')
end
if p_rep < 0.05
     plot(2,mean(pre_rep)+2*std(pre_rep),'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Saccade Amplitude (dva)')
title(sprintf(['ANOVA Saccades 3-5: \n p_{lesion} = ' num2str(P_ANOVA_sacamps(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_sacamps(1),2) ', p_{inter} = ' num2str(P_ANOVA_sacamps(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')


norm_pre_nov_sacamp = pre_sacamp{1}(:,1:median_pre_nov_sac_count);
norm_pre_rep_sacamp = pre_sacamp{2}(:,1:median_pre_rep_sac_count);
norm_pre_rep_sacamp = norm_pre_rep_sacamp./nanmean(norm_pre_rep_sacamp(:,1));
norm_pre_nov_sacamp = norm_pre_nov_sacamp./nanmean(norm_pre_nov_sacamp(:,1));
norm_post_nov_sacamp = post_sacamp{1}(:,1:median_post_nov_sac_count);
norm_post_rep_sacamp = post_sacamp{2}(:,1:median_post_rep_sac_count);
norm_post_rep_sacamp = norm_post_rep_sacamp./nanmean(norm_post_rep_sacamp(:,1));
norm_post_nov_sacamp = norm_post_nov_sacamp./nanmean(norm_post_nov_sacamp(:,1));


norm_pre_nov = mean(norm_pre_nov_sacamp(:,3:15),2);
norm_pre_rep = mean(norm_pre_rep_sacamp(:,3:15),2);
norm_post_nov = mean(norm_post_nov_sacamp(:,3:15),2);
norm_post_rep = mean(norm_post_rep_sacamp(:,3:15),2);

vals = [norm_pre_nov; norm_pre_rep; norm_post_nov; norm_post_rep];
[P_ANOVA_sacamps] = anovan(vals,ids,'model','interaction','varnames',{'Novelty','Lesion'});
%vartestn(vals,ids(:,2))%test if lesion effects variance measures

subplot(2,2,2)
hold all
errorbar(nanmean(norm_pre_nov_sacamp),nanstd(norm_pre_nov_sacamp)./sqrt(npre))
errorbar(nanmean(norm_pre_rep_sacamp),nanstd(norm_pre_rep_sacamp)./sqrt(npre))
errorbar(nanmean(norm_post_nov_sacamp),nanstd(norm_post_nov_sacamp)./sqrt(npost))
errorbar(nanmean(norm_post_rep_sacamp),nanstd(norm_post_rep_sacamp)./sqrt(npost))
hold off
xlabel('Ordinal Saccade #')
ylabel('Relative Saccade Amplitude')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat','Location','NorthEastOutside')
xlim([0 median_post_nov_sac_count+1])
title('Normalized Saccade Amplitudes')


[~,norm_p_nov] = ttest2(norm_pre_nov,norm_post_nov);
[~,norm_p_rep] = ttest2(norm_pre_rep,norm_post_rep);

subplot(2,2,4)
errorb([[mean(norm_pre_nov);mean(norm_post_nov)] [mean(norm_pre_rep);mean(norm_post_rep)]]',...
    [[std(norm_pre_nov)./sqrt(npre) ;std(norm_post_nov)./sqrt(npost)] ...
    [std(norm_pre_rep)./sqrt(npre) ;std(norm_post_rep)./sqrt(npost)]]')
hold on
if norm_p_nov < 0.05
    plot(1,mean(norm_pre_nov)+2*std(norm_pre_nov),'k*')
end
if norm_p_rep < 0.05
     plot(2,mean(norm_pre_rep)+2*std(norm_pre_rep),'k*')
end
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Novel','Repeat'})
ylabel('Relative Saccade Amplitude')
title(sprintf(['ANOVA Saccades 3-5: \n p_{lesion} = ' num2str(P_ANOVA_sacamps(2),2) ...
    ' , p_{Novelty} = ' num2str(P_ANOVA_sacamps(1),2) ', p_{inter} = ' num2str(P_ANOVA_sacamps(3),2)]))
axis square
box off
legend('Pre','Post','Location','NorthEastOutside')

subtitle([pre_files{1}(1:2) ':  Pre vs. Post SCMRM'])

%%
%---Plot Pupil Diameter Over Time---%
t = 1:5:7500;
figure
hold on
for c = 1:2
    dofill(t,pre_pupil{c}/1000,clrs(c),1,60)
end
for c = 1:2
    dofill(t,post_pupil{c}/1000,clrs(c+2),1,60)
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
title('Pupil Diameter')
legend('Pre-Novel','Pre-Repeat','Post-Novel','Post-Repeat')


subtitle([pre_files{1}(1:2) ':  Pre vs. Post SCMRM'])
