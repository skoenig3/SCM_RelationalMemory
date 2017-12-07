%draft code to look at how eye movement patterns as assessed by the recurrence
%analysis predict perofrmance during the SCMRM task pre-lesion.
%written by Seth Konig.
clar

files = {'TO150901.2','TO150903.2','TO150904.3',...
    'TO150908.2','TO150909.2','TO150910.2','TO150911.2',...
    'TO150915.2','TO150916.2','TO150918.2',...
    'TO151001.2','TO151002.2','TO151006.2','TO151007.2',...
    'PW150901.2','PW150902.2','PW150903.2','PW150904.2',...
    'PW150908.2','PW150910.2','PW150911.2',...
    'PW150914.2','PW150916.2','PW150917.2','PW150918.2',...
    'PW150921.2','PW150922.2','PW150923.2','PW150924.2',...
    'MF170222.2','MF170223.2','MF170224.2','MF170227.2',...
    'MF170301.2','MF170303.2','MF170306.2','MF170307.2',...
    'MF170308.2','MF170309.2','MF170313.2','MF170314.2',...
    'MF170315.2','MF170316.2','MF170317.2',...
    'RR160425.2','RR160426.2','RR160701.2','RR160706.2','RR160707.2',...
    'RR160708.2','RR160714.2','RR160715.2','RR160718.2','RR160721.2',...
    'RR160725.2','RR160726.2','RR160727.2','RR160728.2','RR160729.2',...
    'RR160801.2','RR160803.2','RR160804.2','RR160805.2'};

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
all_corm = NaN(1,length(files)*36);
mean_fix_durs = NaN(1,length(files)*36);
mean_sac_amps = NaN(1,length(files)*36);

all_nov_fix_durs = NaN(length(files)*36,40);
all_rep_fix_durs = NaN(length(files)*36,40);
all_nov_sac_amps = NaN(length(files)*36,40);
all_rep_sac_amps = NaN(length(files)*36,40);
all_set_nums =  NaN(1,length(files)*36);
which_monkey = NaN(1,length(files)*36);
all_ROI_time = NaN(1,length(files)*36);

all_nov_recurrence_maps = cell(1,2300);
for n = 1:2300;
    all_nov_recurrence_maps{n} = NaN(N_fix,N_fix);
end

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
        
        if nov_img_off-nov_img_on > 10500
            continue
        end
        
        nov_x = fixationstats{2*img-1}.XY(1,nov_img_on:nov_img_off);
        nov_y = fixationstats{2*img-1}.XY(2,nov_img_on:nov_img_off);
        
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
        nov_sactimes = nov_sactimes-nov_img_on;
        
        
        
        %for repeat images
        rep_allval = per(2*img).allval;
        rep_alltim = per(2*img).alltim;
        rep_img_on = rep_alltim(rep_allval == 23)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        rep_img_off = rep_alltim(rep_allval == 24)-rep_alltim(rep_allval == 100);%image on relative to eye data start
        
        if rep_img_off-rep_img_on > 10500
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
        rep_sactimes = rep_sactimes-rep_img_on;
        
        
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
%         
%         if rep_trialtype ==4 %moved image. Data from the old ROI (where the object was before being moved)
%             ROIfix2_oldROI = find(...
%                 (rep_fix(1,:) > ROI1(1) & rep_fix(1,:) < ROI1(2) & ...
%                 rep_fix(2,:) > ROI1(3) & rep_fix(2,:) < ROI1(4)));%ROIfix2_oldROI contains indices of rep_fix corresponding to fixations that occur in the 'old' ROI of the 'moved' image
%         else
%             ROIfix2_oldROI = [];
%         end
        
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
%             for j=1:length(ROIfix2_oldROI)
%                 tempvec3(rep_fixtimes(1,ROIfix2_oldROI(j)):rep_fixtimes(2,ROIfix2_oldROI(j)))=1;
%             end
%             
            %---Calculate Stats for the Novel image---%
            nov_propall_fixations = length(ROIfix1)/length(nov_fix);
            nov_prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
            nov_propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
            nov_prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
            nov_prop15_time = nansum(tempvec1(500:1500))/sum(~isnan(tempvec1(500:1500)));
            
            %---Calculate Stats for the Second Presentation---%
            if trialtype == 4 %for moved
                rep_propall_time = (nansum(tempvec2)+nansum(tempvec3))/sum(~isnan(tempvec2));
            else %for replaced or repeat
                rep_propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
            end
        else
            continue
        end
        
        
        %%
        %fixation duration by ordinal fixation number
        fixdurs = nov_fixtimes(2,:)-nov_fixtimes(1,:)+1;
        repfixdurs = rep_fixtimes(2,:)-rep_fixtimes(1,:)+1;
        
        [recurrence_rate,recurrence_map,corm,laminarity,laminarity_len,...
            forward_trace,forward_trace_len,reverse_trace,reverse_trace_len] = ...
            calculate_auto_recurrence(nov_fix,N_fix,dist_threshold);
        
        img_id = img_id+1;
        all_set_nums(img_id) = sets;
        all_img_types(img_id) = rep_trialtype;
        all_recurrence_rates(img_id) = recurrence_rate;
        all_laminarity(img_id) = laminarity;
        all_forward_trace(img_id) = forward_trace;
        all_reverse_trace(img_id) =reverse_trace;
        all_corm(img_id) = corm;
        
        all_ROI_time(img_id) = rep_propall_time;
        mean_fix_durs(img_id) = mean(fixdurs(4:end));
        
        all_nov_recurrence_maps{img_id} = recurrence_map;
        
        if length(fixdurs) < 40
            all_nov_fix_durs(img_id,1:length(fixdurs)) = fixdurs;
        else
            all_nov_fix_durs(img_id,:) = fixdurs(1:40);
        end
        
        if length(repfixdurs) < 40
            all_rep_fix_durs(img_id,1:length(repfixdurs)) = repfixdurs;
        else
            all_rep_fix_durs(img_id,:) = repfixdurs(1:40);
        end
        
        sac_amps = NaN(1,size(nov_sactimes,2));
        for s = 1:size(nov_sactimes,2)
            sac_amps(s) = sqrt((nov_x(nov_sactimes(2,s))-nov_x(nov_sactimes(1,s))).^2 + ...
                (nov_y(nov_sactimes(2,s))-nov_y(nov_sactimes(1,s))).^2)/24;
        end
        if length(sac_amps) < 40
            all_nov_sac_amps(img_id,1:length(sac_amps)) = sac_amps;
        else
            all_nov_sac_amps(img_id,:) = sac_amps(1:40);
        end
        
        sac_amps = NaN(1,size(rep_sactimes,2));
        for s = 1:size(rep_sactimes,2)
            sac_amps(s) = sqrt((rep_x(rep_sactimes(2,s))-rep_x(rep_sactimes(1,s))).^2 + ...
                (rep_y(rep_sactimes(2,s))-rep_y(rep_sactimes(1,s))).^2)/24;
        end
        if length(sac_amps) < 40
            all_rep_sac_amps(img_id,1:length(sac_amps)) = sac_amps;
        else
            all_rep_sac_amps(img_id,:) = sac_amps(1:40);
        end
        
        mean_sac_amps(img_id) = nanmean(sac_amps(4:end));
        
        if strcmp('TO',files{sets}(1:2))
            which_monkey(img_id) = 1;
        elseif strcmp('PW',files{sets}(1:2))
            which_monkey(img_id) = 2;
        elseif strcmp('MF',files{sets}(1:2))
            which_monkey(img_id) = 3;
        elseif strcmp('RR',files{sets}(1:2))
            which_monkey(img_id) = 4;
        end
        
    end
end
%%
%normalize fixation durations to first fixation
for monk = 1:4
    these_trials = find(which_monkey == monk);
    nov_fix_dur = nanmean(all_nov_fix_durs(these_trials,1));
    all_nov_fix_durs(these_trials,:) = all_nov_fix_durs(these_trials,:)./nov_fix_dur;
    all_rep_fix_durs(these_trials,:) = all_rep_fix_durs(these_trials,:)./nov_fix_dur;
end
%%

good_nov_durs = [];
good_rep_durs = [];
bad_nov_durs = [];
bad_rep_durs = [];

good_nov_amps = [];
good_rep_amps = [];
bad_nov_amps = [];
bad_rep_amps = [];

all_dprimes = cell(1,3);
combos = triu(ones(7,7));

figure
for type = 2:4
    all_dprimes{type-1} = zeros(3,7);
    
    these_recurrence_rates = all_recurrence_rates(all_img_types == type);
    these_laminarity = all_laminarity(all_img_types == type);
    these_fwd_retraces = all_forward_trace(all_img_types == type);
    these_rev_retraces = all_reverse_trace(all_img_types == type);
    these_corm = all_corm(all_img_types == type);
    these_fix_durs = mean_fix_durs(all_img_types == type);
    these_sac_amps = mean_sac_amps(all_img_types == type);
    
    var_names = {'these_recurrence_rates','these_laminarity','these_fwd_retraces','these_rev_retraces',...
        'these_corm','these_fix_durs','these_sac_amps'};
 
    %have to do columns and rows seperately to get all combos and itself
    %seperately
    for n1 = 1:7
        use_these_vars = var_names(find(combos(n1,:)));
        group_vals = [];
        for u = 1:length(use_these_vars)
            group_vals = [group_vals eval(use_these_vars{u})'];
        end
        [U,~,~] = pca(group_vals,1);
        pct33 = prctile(U,33);
        pct67 = prctile(U,67);
        bad_ROI = group_ROI(U <= pct33);
        good_ROI = group_ROI(U >= pct67);
        all_dprimes{type-1}(1,n1) = abs(mean(good_ROI)-mean(bad_ROI))/((std(good_ROI)+std(bad_ROI))/2);
        
        
        use_these_vars = var_names(find(combos(:,n1)));
        group_vals = [];
        for u = 1:length(use_these_vars)
            group_vals = [group_vals eval(use_these_vars{u})'];
        end
        [U,~,~] = pca(group_vals,1);
        pct33 = prctile(U,33);
        pct67 = prctile(U,67);
        bad_ROI = group_ROI(U <= pct33);
        good_ROI = group_ROI(U >= pct67);
        all_dprimes{type-1}(2,n1) = abs(mean(good_ROI)-mean(bad_ROI))/((std(good_ROI)+std(bad_ROI))/2);
        
        use_these_vars = var_names(n1);
        U = eval(use_these_vars{1});
        pct33 = prctile(U,33);
        pct67 = prctile(U,67);
        bad_ROI = group_ROI(U <= pct33);
        good_ROI = group_ROI(U >= pct67);
        all_dprimes{type-1}(3,n1) = abs(mean(good_ROI)-mean(bad_ROI))/((std(good_ROI)+std(bad_ROI))/2);
    end
    
    
    %group_vals = [all_recurrence_rates(all_img_types == type)' all_corm(all_img_types == type)'...
 %       all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...
%         all_reverse_trace(all_img_types == type)' mean_fix_durs(all_img_types == type)'];%...
    %mean_sac_amps(all_img_types == type)'];
    
%      group_vals = [all_corm(all_img_types == type)' mean_fix_durs(all_img_types == type)' ...
%        all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...

    %all_laminarity(all_img_types == type)'...
    %all_recurrence_rates(all_img_types == type)'
    %     group_vals = [all_recurrence_rates(all_img_types == type)' all_corm(all_img_types == type)'...
    %         all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...
    %         all_reverse_trace(all_img_types == type)' mean_fix_durs(all_img_types == type)'];
    group_ROI = all_ROI_time(all_img_types == type);
    all_nov_group_fix_durs = all_nov_fix_durs(all_img_types == type,:);
    all_rep_group_fix_durs = all_rep_fix_durs(all_img_types == type,:);
    all_nov_group_sac_amps = all_nov_sac_amps(all_img_types == type,:);
    all_rep_group_sac_amps = all_rep_sac_amps(all_img_types == type,:);
    
    these_maps = all_nov_recurrence_maps(all_img_types == type);
    
    [U,~,~] = pca(group_vals,1);
    pct33 = prctile(U,33);
    pct67 = prctile(U,67);
    
    bad_ROI = group_ROI(U <= pct33);
    good_ROI = group_ROI(U >= pct67);
    
     good_fwd_retraces = these_fwd_retraces(U >= pct67);
    bad_fwd_retraces = these_fwd_retraces(U <= pct33);
    
       good_rev_retraces = these_rev_retraces(U >= pct67);
    bad_rev_retraces = these_rev_retraces(U <= pct33);
    
    
    these_corm = all_corm(all_img_types == type);
    good_corm = these_corm(U >= pct67);
    bad_corm = these_corm(U <= pct33);

    good_rates = these_recurrence_rates(U >= pct67);
    bad_rates = these_recurrence_rates(U <= pct33);
    
    good_maps = these_maps(U >= pct67);
    good_maps = nanmean(reshape(cell2mat(good_maps),[40,40,length(good_maps)]),3);
    bad_maps = these_maps(U <= pct33);
    bad_maps = nanmean(reshape(cell2mat(bad_maps),[40,40,length(bad_maps)]),3);
    
    figure(10+type)
    subplot(2,2,1)
    gm = plot_average_recurrence_plots(good_maps,0,25);
    title('Good')
    
    subplot(2,2,2)
    bm = plot_average_recurrence_plots(bad_maps,0,25);
    title('Bad')
    
    subplot(2,2,3)
    plot_difference_of_average_recurrence_plots(bm,gm)
    title('Good-Bad')
    
    if type == 2
        subtitle(['Repeat, p_{whole} = ' num2str(p,2)])
    elseif type == 3
        subtitle(['Replaced, p_{whole} = ' num2str(p,2)])
    else
        subtitle(['Moved, p_{whole} = ' num2str(p,2)])
    end

    good_nov_durs = [good_nov_durs; all_nov_group_fix_durs(U >= pct67,:)];
    good_rep_durs = [good_rep_durs; all_rep_group_fix_durs(U >= pct67,:)];
    bad_nov_durs = [bad_nov_durs;   all_nov_group_fix_durs(U <= pct33,:)];
    bad_rep_durs = [bad_rep_durs;   all_rep_group_fix_durs(U <= pct33,:)];
    
    good_nov_amps = [good_nov_amps; all_nov_group_sac_amps(U >= pct67,:)];
    good_rep_amps = [good_rep_amps; all_rep_group_sac_amps(U >= pct67,:)];
    bad_nov_amps = [bad_nov_amps;   all_nov_group_sac_amps(U <= pct33,:)];
    bad_rep_amps = [bad_rep_amps;   all_rep_group_sac_amps(U <= pct33,:)];
    
    [~,p] = ttest2(good_ROI,bad_ROI);
    
    figure(101)
    subplot(2,2,type-1)
    hold on
    bar([mean(bad_ROI) mean(good_ROI)])
    errorb([mean(bad_ROI) mean(good_ROI)],[std(bad_ROI)/sqrt(length(bad_ROI)) std(good_ROI)/sqrt(length(good_ROI))])
    set(gca,'Xtick',[1:2])
    set(gca,'XtickLabel',{'Good','Bad'})
    ylabel('% of Time in ROI')
    axis square
    
    if type == 2
        title(['Repeat, p_{whole} = ' num2str(p,2)])
    elseif type == 3
        title(['Replaced, p_{whole} = ' num2str(p,2)])
    else
        title(['Moved, p_{whole} = ' num2str(p,2)])
    end
end

subplot(2,2,4)
hist(U,25)
xlabel('PCA Values')
ylabel('Count')
axis square
box off

subtitle('PCA sorted')

type = 2;
group_vals = [all_recurrence_rates(all_img_types == type)' all_corm(all_img_types == type)'...
    all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...
    all_reverse_trace(all_img_types == type)' mean_fix_durs(all_img_types == type)'];% ...
%mean_sac_amps(all_img_types == type)'];
%     group_vals = [all_corm(all_img_types == type)' mean_fix_durs(all_img_types == type)' ...
%         all_forward_trace(all_img_types == type)' all_laminarity(all_img_types == type)'...
%         all_reverse_trace(all_img_types == type)'];
group_ROI = all_ROI_time(all_img_types == type);

[U,S,V] = pca(group_vals,1);
pct33 = prctile(U,40);
pct67 = prctile(U,60);

bad_ROI = group_ROI(U <= pct33);
good_ROI = group_ROI(U >= pct67);

N=sampsizepwr('t',[mean(good_ROI) std(good_ROI)],[mean(bad_ROI) std(bad_ROI)],0.8)

figure
subplot(2,2,1)
hold on
plot(nanmean(all_nov_fix_durs(:,1:20)))
plot(nanmean(all_rep_fix_durs(:,1:20)))
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('All Image Presentations')
box off

subplot(2,2,2)
hold on
plot(nanmean(good_nov_durs(:,1:20)))
plot(nanmean(good_rep_durs(:,1:20)))
plot(nanmean(bad_nov_durs(:,1:20)))
plot(nanmean(bad_rep_durs(:,1:20)))
hold off
xlabel('Ordinal Fixation #')
ylabel('Fixation Duration (ms)')
title('All Image Presentations')
box off
legend('Hi Nov','Hi Rep','Low Nov','Low Rep')

good_nov_vals = nanmean(good_nov_durs,2);
good_rep_vals = nanmean(good_rep_durs,2);
bad_nov_vals = nanmean(bad_nov_durs,2);
bad_rep_vals = nanmean(bad_rep_durs,2);

subplot(2,2,3)
hold on
bar([nanmean(good_nov_vals) nanmean(good_rep_vals) nanmean(bad_nov_vals) nanmean(bad_rep_vals)])
errorb([nanmean(good_nov_vals) nanmean(good_rep_vals) nanmean(bad_nov_vals) nanmean(bad_rep_vals)],...
    [nanstd(good_nov_vals) nanstd(good_rep_vals) nanstd(bad_nov_vals) nanstd(bad_rep_vals)]...
    ./sqrt(length(15)))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Hi Nov','Hi Rep','Low Nov','Low Rep'})
box off
xlabel('Presentation Type')
ylabel('Fixation Duration (ms)')
ylim([0.85 1.6])

subtitle('Effect of Fixation Durations on Memory in SCM')


figure
subplot(2,2,1)
hold on
plot(nanmean(all_nov_sac_amps(:,1:20)))
plot(nanmean(all_rep_sac_amps(:,1:20)))
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
title('All Image Presentations')
box off

subplot(2,2,2)
hold on
plot(nanmean(good_nov_amps(:,1:20)))
plot(nanmean(good_rep_amps(:,1:20)))
plot(nanmean(bad_nov_amps(:,1:20)))
plot(nanmean(bad_rep_amps(:,1:20)))
hold off
xlabel('Ordinal Saccade #')
ylabel('Saccade Amplitude (dva)')
title('All Image Presentations')
box off
legend('Hi Nov','Hi Rep','Low Nov','Low Rep')

good_nov_vals = nanmean(good_nov_amps,2);
good_rep_vals = nanmean(good_rep_amps,2);
bad_nov_vals = nanmean(bad_nov_amps,2);
bad_rep_vals = nanmean(bad_rep_amps,2);

subplot(2,2,3)
hold on
bar([nanmean(good_nov_vals) nanmean(good_rep_vals) nanmean(bad_nov_vals) nanmean(bad_rep_vals)])
errorb([nanmean(good_nov_vals) nanmean(good_rep_vals) nanmean(bad_nov_vals) nanmean(bad_rep_vals)],...
    [nanstd(good_nov_vals) nanstd(good_rep_vals) nanstd(bad_nov_vals) nanstd(bad_rep_vals)]...
    ./sqrt(length(15)))
hold off
set(gca,'Xtick',1:4)
set(gca,'XtickLabel',{'Hi Nov','Hi Rep','Low Nov','Low Rep'})
box off
xlabel('Presentation Type')
ylabel('Saccade Amplitude')

subtitle('Effect of Saccade Amplitude on Memory in SCM')

%%
figure