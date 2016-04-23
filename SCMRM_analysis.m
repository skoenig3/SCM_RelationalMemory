% analyze pilot data from Wilbur running various forms of the SCene
% manipulation task
% originally written 3/1/15 Seth Konig
% based somewhat on  "SCM_ROI_FixationTimeOnly_V2_perSession.m".
%%
clear,clc

sets = 1;
cortex_files = {'PW150512.2'};
imgdur = 7;

eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\';


figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';

% for file = 1:length(cortex_files);
%     get_SCM_RM_eyedat(cortex_files{file},sets(file))
% end

imageX = 800;
imageY = 600;

allnumFixationsInROI = cell(length(cortex_files),4); %number of fixations in ROI by session by type
num10FixationsInROI = cell(length(cortex_files),4); %number of fixations in ROI by session by type for 1st 10 fixations

allTime_ROI_timeWindow =  cell(length(cortex_files),4);%amount of time in ROI by session by type
Time3_ROI_timeWindow =  cell(length(cortex_files),4);%amount of time in ROI by session by type in 1st 3 secs

area = cell(length(cortex_files),4);

for SET = 1:length(cortex_files)
    setnum = sets(SET);
    scm_image_dir = ['R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\Image Sets\'...
        'Set0' num2str(sets(SET)) '\'];
    
    load([ROI_dir 'SCM_RM_ROIs_Set' num2str(setnum) '.mat']);%ROIs for image set
    load([eye_data_dir cortex_files{SET}(1:8) '_' cortex_files{SET}(end) '-fixation.mat']);%preprocessed cortex data
    
    all_time = zeros(5,7000);%novel, familiar, replaced, moved,new_moved
    all_time_points = zeros(5,7000);
    for img = 1:36;
        img_ind = find(which_image == img);
        nov_ind = img_ind(1);%1st presentation (nov)
        rep_ind = img_ind(2);%2nd presentation  (rep)
        
        nov_time_out = zeros(1,7000); %time outside image or crosshair fixation
        novx = 24*fixationstats{nov_ind}.XY(1,1:7000/5)+400;
        novy = 24*fixationstats{nov_ind}.XY(2,1:7000/5)+300;
        novfixations = fixationstats{nov_ind}.fixations;
        novfixations(1,:) =  24*novfixations(1,:)+400;
        novfixations(2,:) =  (24*novfixations(2,:)+300);
        novtimes =  fixationstats{nov_ind}.fixationtimes;
        novfixations(2,:) = novfixations(2,:);
        if (novfixations(1,1) > 300 && novfixations(1,1) < 500) && ...
                (novfixations(2,1) > 200 && novfixations(2,1) < 400)
            nov_time_out(1:novtimes(2,1)) = 1;
            novx(1:novtimes(2,1)) = [];
            novy(1:novtimes(2,1)) = [];
            novfixations(:,1) = [];
            novtimes(:,1) = [];
        end
        nov_outside = find(novx > imageX | novx < 1 | novy > imageX | novy < 1);
        nov_time_out(nov_outside) = 1;
        novtimes = novtimes*5; 
        
        rep_time_out = zeros(1,7000);
        repfixations = fixationstats{rep_ind}.fixations;
        repfixations(1,:) =  24*repfixations(1,:)+400;
        repfixations(2,:) = (24*repfixations(2,:)+300);
        reptimes =  fixationstats{rep_ind}.fixationtimes;
        repfixations(2,:) = repfixations(2,:);
        repx = 24*fixationstats{rep_ind}.XY(1,1:7000/5)+400;
        repy = 24*fixationstats{rep_ind}.XY(2,1:7000/5)+300;
        if (repfixations(1,1) > 300 && repfixations(1,1) < 500) && ...
                (repfixations(2,1) > 200 && repfixations(2,1) < 400)
            rep_time_out(1:reptimes(2,1)) = 1;
            repfixations(:,1) = [];
            reptimes(:,1) = [];
            repx(1:reptimes(2,1)) = [];
            repy(1:reptimes(2,1)) = [];
        end
        rep_outside = find(repx > imageX | repx < 1 | repy > imageX | repy < 1);
        rep_time_out(rep_outside) = 1;
        reptimes = reptimes*5; 
        
        novy(novx > imageX) = [];
        novx(novx > imageX) = [];
        novx(novy > imageY) = [];
        novy(novy > imageY) = [];
        novy(novx < 1) = [];
        novx(novx < 1) = [];
        novx(novy < 1) = [];
        novy(novy < 1) = [];
        
        repy(repx > imageX) = [];
        repx(repx > imageX) = [];
        repx(repy > imageY) = [];
        repy(repy > imageY) = [];
        repy(repx < 1) = [];
        repx(repx < 1) = [];
        repx(repy < 1) = [];
        repy(repy < 1) = [];
        
        %ignore anything after 7000 ms
        late = find(novtimes(1,:) > 7000);
        novtimes(:,late) = [];
        novfixations(:,late) = [];
        novtimes(novtimes > 7000) = 7000;
        
        late = find(reptimes(1,:) > 7000);
        reptimes(:,late) = [];
        repfixations(:,late) = [];
        reptimes(reptimes > 7000) = 7000;
        
        if trialtype(rep_ind) == 2 || trialtype(rep_ind) == 3 %familiar and replaced, respectively
            ROI1 = ROIs{img}{1};
            ROI1(1) = ROI1(1)-24;
            ROI1(2) = ROI1(2)+24;
            ROI1(3) = ROI1(3)-24;
            ROI1(4) = ROI1(4)+24;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 =ROI1;
        elseif trialtype(rep_ind) == 4 %moved
            ROI1 = ROIs{img}{1};%original location
            ROI1(1) = ROI1(1)-24;
            ROI1(2) = ROI1(2)+24;
            ROI1(3) = ROI1(3)-24;
            ROI1(4) = ROI1(4)+24;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 = ROIs{img}{2};%new location after move
            ROI2(1) = ROI2(1)-24;
            ROI2(2) = ROI2(2)+24;
            ROI2(3) = ROI2(3)-24;
            ROI2(4) = ROI2(4)+24;
            ROI2(ROI2 < 1) = 1;
            ROI2(ROI2 > imageX) = imageX;
            if ROI2(3) > imageY; ROI2(3) = imageY;end
            if ROI2(4) > imageY; ROI2(4) = imageY;end
        end
        
        %added commented out code below to ploty eye traces on
        %images.
        
        
        rep_trialType = trialtype(rep_ind);
        
        ROIfix1 = find(...
            (novfixations(1,:) > ROI1(1) & novfixations(1,:) < ROI1(2)) & ...
            (novfixations(2,:) > ROI1(3) & novfixations(2,:) < ROI1(4))); %ROIfix1 contains indices of novfixations corresponding to fixations that occur in ROI of novel image.
        ROIfix2 = find(...
            repfixations(1,:) > ROI2(1) & repfixations(1,:) < ROI2(2) & ...
            repfixations(2,:) > ROI2(3) & repfixations(2,:) < ROI2(4));%ROIfix2 contains indices of repfixations corresponding to fixations that occur in ROI of the manipulated image
  
        if trialtype(rep_ind) ==4 %moved image. Data from the old ROI (where the object was before being moved)
            ROIfix2_oldROI = find(...
                (repfixations(1,:) > ROI1(1) & repfixations(1,:) < ROI1(2) & ...
                repfixations(2,:) > ROI1(3) & repfixations(2,:) < ROI1(4)));%ROIfix2_oldROI contains indices of repfixations corresponding to fixations that occur in the 'old' ROI of the 'moved' image
        else
            ROIfix2_oldROI = [];
        end
        if ~isempty(ROIfix2_oldROI)
            disp('found fixations')
        end
%         figure
%         subplot(2,2,1)
%         hold on
%         image(flipud(imread([scm_image_dir num2str(which_image(nov_ind)) '.bmp'])));
%         plot(novx,novy);
%         plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
%             [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
%         plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
%             [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
%         hold off
%         box off
%         xlim([0 imageX])
%         ylim([0 imageY])
%         axis off
%         axis equal
%         
%         subplot(2,2,3)
%         hold on
%         if rep_trialType == 2 % familiar
%             image(flipud(imread([scm_image_dir num2str(which_image(rep_ind)) 'p.bmp'])));
%         elseif rep_trialType ==3 %replaced
%             image(flipud(imread([scm_image_dir num2str(which_image(rep_ind)) 'r.bmp'])));
%         elseif rep_trialType ==4 %moved
%             image(flipud(imread([scm_image_dir num2str(which_image(rep_ind)) 'm.bmp'])));
%         end
%         plot(repx,repy);
%         plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
%             [ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
%         plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
%             [ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
%         hold off
%         box off
%         xlim([0 imageX])
%         ylim([0 imageY])
%         axis off
%         axis equal
        
         if isempty(ROIfix1)%didnt look at original region, moving on to next image
%             subplot(2,2,2)
%             title('No fixations in Novel ROI')
%             save_and_close_fig(figure_dir,num2str(img));
            continue
        end
        
        
        
        %%%%Code for calculating avg time in ROI for each image type
        tempvec1=zeros(1,7000);%novel image
        for j=1:length(ROIfix1)
            tempvec1(novtimes(1,ROIfix1(j)):novtimes(2,ROIfix1(j)))=1;
        end
        tempvec1(find(nov_time_out)) = NaN;
        
        tempvec2=zeros(1,7000);%manipulated image new ROI
        for j=1:length(ROIfix2)
            tempvec2(reptimes(1,ROIfix2(j)):reptimes(2,ROIfix2(j)))=1;
        end
        tempvec2(find(rep_time_out)) = NaN;
        
        tempvec3=zeros(1,7000);%manipulated imaged old moved ROI
        for j=1:length(ROIfix2_oldROI)
            tempvec3(reptimes(1,ROIfix2_oldROI(j)):reptimes(2,ROIfix2_oldROI(j)))=1;
        end
        tempvec3(find(rep_time_out)) = NaN;
        
        area_ROI1 = abs(ROI1(1)-ROI1(2))*abs(ROI1(3)-ROI1(4)); %area of ROI for novel image
        area_ROI2 = abs(ROI2(1)-ROI2(2))*abs(ROI2(3)-ROI2(4)); %for replaced, moved or familiar
        
        %for novel presentation
        area{SET,1} = [area{SET,1} area_ROI1];
        
        propall_fixations = length(ROIfix1)/length(novfixations);
        prop10_fixaitons = length(ROIfix1(ROIfix1 <= 10))/10;
        propall_time = nansum(tempvec1)./sum(~isnan(tempvec1));
        prop3_time = nansum(tempvec1(1:3000))/sum(~isnan(tempvec1(1:3000)));
        
        allnumFixationsInROI{SET,1} =[allnumFixationsInROI{SET,1} propall_fixations];
        num10FixationsInROI{SET,1} = [num10FixationsInROI{SET,1} prop10_fixaitons];
        allTime_ROI_timeWindow{SET,1} = [allTime_ROI_timeWindow{SET,1} propall_time];
        Time3_ROI_timeWindow{SET,1} = [Time3_ROI_timeWindow{SET,1} prop3_time];
        
%         s(1) = subplot(2,2,2);
%         bar([propall_fixations propall_time prop10_fixaitons prop3_time])
%         set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
%         ylabel('Proportion in ROI')
%         title(['First Fixation in ROI is # ' num2str(ROIfix1(1))])
        
        
        %for second presentation
        if trialtype(rep_ind) == 4
            area{SET,trialtype(rep_ind)} = [area{SET,trialtype(rep_ind)} area_ROI2+area_ROI1];
            
               
            propall_fixations = (length(ROIfix2)+length(ROIfix2_oldROI))/length(repfixations);
            prop10_fixaitons = (length(ROIfix2(ROIfix2 <= 10))+...
              length(ROIfix2_oldROI(ROIfix2_oldROI <= 10)))./10;
            propall_time = (nansum(tempvec2)+nansum(tempvec3))/sum(~isnan(tempvec2));
            prop3_time = (nansum(tempvec2(1:3000))+nansum(tempvec3(1:3000)))...
                /sum(~isnan(tempvec2(1:3000)));
            
            allnumFixationsInROI{SET,trialtype(rep_ind)} =[allnumFixationsInROI{SET,trialtype(rep_ind)} propall_fixations];
            num10FixationsInROI{SET,trialtype(rep_ind)} = [num10FixationsInROI{SET,trialtype(rep_ind)} prop10_fixaitons];
            allTime_ROI_timeWindow{SET,trialtype(rep_ind)} = [allTime_ROI_timeWindow{SET,trialtype(rep_ind)} propall_time];
            Time3_ROI_timeWindow{SET,trialtype(rep_ind)} = [Time3_ROI_timeWindow{SET,trialtype(rep_ind)} prop3_time];
            
%             s(2) = subplot(2,2,4);
%             bar([propall_fixations propall_time prop10_fixaitons prop3_time])
%             set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
%             ylabel('Proportion in ROI')
%             if ~isempty(ROIfix2)
%                 title(['First Fixation in ROI is # ' num2str(ROIfix2(1))])
%             end
            
        else %for replaced or repeat
            area{SET,trialtype(rep_ind)} = [area{SET,trialtype(rep_ind)} area_ROI2];
            
            propall_fixations = length(ROIfix2)/length(repfixations);
            prop10_fixaitons = length(ROIfix2(ROIfix2 <= 10))/10;
            propall_time = nansum(tempvec2)./sum(~isnan(tempvec2));
            prop3_time = nansum(tempvec2(1:3000))/sum(~isnan(tempvec2(1:3000)));
            
            allnumFixationsInROI{SET,trialtype(rep_ind)} =[allnumFixationsInROI{SET,trialtype(rep_ind)} propall_fixations];
            num10FixationsInROI{SET,trialtype(rep_ind)} = [num10FixationsInROI{SET,trialtype(rep_ind)} prop10_fixaitons];
            allTime_ROI_timeWindow{SET,trialtype(rep_ind)} = [allTime_ROI_timeWindow{SET,trialtype(rep_ind)} propall_time];
            Time3_ROI_timeWindow{SET,trialtype(rep_ind)} = [Time3_ROI_timeWindow{SET,trialtype(rep_ind)} prop3_time];
            
%             s(2) = subplot(2,2,4);
%             bar([propall_fixations propall_time prop10_fixaitons prop3_time])
%             set(gca,'XtickLabel',{'Fixations_{all}','Time_{all}','Fixations_{10}','Time_3'});
%             ylabel('Proportion in ROI')
%             if ~isempty(ROIfix2)
%                 title(['First Fixation in ROI is # ' num2str(ROIfix2(1))])
%             end
            
        end
%         yl1 = get(s(1),'ylim');
%         yl2 = get(s(2),'ylim');
%         ymax = max([yl1 yl2]);
%         set(s(1),'ylim',[0 ymax]);
%         set(s(2),'ylim',[0 ymax]);
%         _and_close_fig(figure_dir,num2str(img));
        
        
        all_time_points(1,~isnan(tempvec1)) = all_time_points(1,~isnan(tempvec1))+1;
        tempvec1(isnan(tempvec1)) = 0;
        all_time(1,:) = all_time(1,:) + tempvec1;
        
        all_time_points(trialtype(rep_ind),~isnan(tempvec2)) = ...
            all_time_points(trialtype(rep_ind),~isnan(tempvec2))+1;
        tempvec2(isnan(tempvec2)) = 0;
        all_time(trialtype(rep_ind),:) =  all_time(trialtype(rep_ind),:)+tempvec2;
        
        if trialtype(rep_ind) == 4 %for moved only
            all_time_points(5,~isnan(tempvec3)) = ...
                all_time_points(5,~isnan(tempvec3))+1;
            tempvec3(isnan(tempvec3)) = 0;
            all_time(5,:) =  all_time(5,:)+tempvec3+tempvec2;
        end
    end
    
    all_time_points(all_time_points == 0) = 1; %for division purposes
    all_time = all_time./all_time_points; %normalize
    
    clrs = 'brmgg';
    line = {'-','-','-','-','--'};
    image_count = cellfun(@length,allnumFixationsInROI(SET,:));%num images in each set
    image_count = [image_count image_count(end)];
    figure
    hold on
    for r = 1:size(all_time,1)
        plot(100*filtfilt(1/200*ones(1,200),1,all_time(r,:)),...
            [clrs(r) line{r}]);
    end
    hold off
    xlabel('Time (ms)')
    ylabel('% of time in ROI')
    legend('Novel','Repeat','Replaced','Moved (new)','Moved (new+old)');
    title(['Set ' num2str(sets(SET))])
end
%%
SET =1 ;
figure
bar(cellfun(@mean,allnumFixationsInROI(SET,:)))
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of all Fixations in ROI')
title(['Set ' num2str(sets(SET))])
%%
figure
bar(cellfun(@mean,num10FixationsInROI(SET,:)))
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of 1st 10 Fixations in ROI')
title(['Set ' num2str(sets(SET))])
%%
figure
bar(cellfun(@mean,num10FixationsInROI(SET,:)))
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of 1st 10 Fixations in ROI')
title(['Set ' num2str(sets(SET))])
%%
figure
bar(cellfun(@mean,allTime_ROI_timeWindow(SET,:)))
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of all Time in ROI')
title(['Set ' num2str(sets(SET))])
%%
figure
bar(cellfun(@mean,Time3_ROI_timeWindow(SET,:)))
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of 1st 3 seconds in ROI')
title(['Set ' num2str(sets(SET))])
