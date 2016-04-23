% analyze pilot data from Wilbur running various forms of the SCene
% manipulation task
% originally written 3/1/15 Seth Konig
% based somewhat on  "SCM_ROI_FixationTimeOnly_V2_perSession.m".

clear,clc

eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\';

sets = [1 2];
cortex_files = {'WR150226.3','WR150227.2'};
imgdur = [5 5];
img_seperation = [12 0];
nmin = 1;
nmax = 10;

% for file = 1:length(cortex_files);
%     get_SCM_RM_eyedat(cortex_files{file},sets(file))
% end

screen_size = get(0, 'ScreenSize');
imageX = 800; imageY = 600;

numFixationsInROI = cell(length(cortex_files),4); %number of fixations in ROI by session by typ
Time_ROI_timeWindow =  cell(length(cortex_files),4);%amount of time in ROI by session by type
prop_fixROI_timeWindow = cell(length(cortex_files),4); %number of fixations in ROI during time interval by session by type
area = cell(length(cortex_files),4);

for SET = 1:length(cortex_files)
    SETNUM = sets(SET);
    
    load([ROI_dir 'Set' num2str(SET) '_ROIs.mat']);%ROIs for image set
    load([eye_data_dir cortex_files{SET}(1:8) '_' cortex_files{SET}(end) '-fixation.mat']);%preprocessed cortex data
    
    all_time = zeros(5,5000);%novel, familiar, replaced, moved,new_moved
    for img = 1:36;
        img_ind = find(which_image == img);
        nov_ind = img_ind(1);%1st presentation (nov)
        rep_ind = img_ind(2);%2nd presentation  (rep)
        novfixations = fixationstats{nov_ind}.fixations;
        novfixations(1,:) =  24*novfixations(1,:)+400;
        novfixations(2,:) =  600-(24*novfixations(2,:)+300);
        novtimes =  5*fixationstats{nov_ind}.fixationtimes;
        novfixations(2,:) = novfixations(2,:);
        repfixations = fixationstats{rep_ind}.fixations;
        repfixations(1,:) =  24*repfixations(1,:)+400;
        repfixations(2,:) = 600-(24*repfixations(2,:)+300);
        reptimes =  5*fixationstats{rep_ind}.fixationtimes;
        repfixations(2,:) = repfixations(2,:);
        
        
        %ignore anything after 5000 ms
        novtimes(novtimes > 5000) = 5000;
        reptimes(reptimes > 5000) = 5000;
        %Calculate the proportion of the first 10 fixations spent in ROI and in per pixel of ROI
        if length(novfixations) > nmax
            novfixations = novfixations(:,nmin:nmax);
        end
        if length(repfixations)> nmax
            repfixations= repfixations(:,nmin:nmax);
        end
        if trialtype(rep_ind) == 2 || trialtype(rep_ind) == 3 %familiar and replaced, respectively
            ROI1 = ROIs{img};
            ROI1(1) = ROI1(1)-24;
            ROI1(2) = ROI1(2)+24;
            ROI1(3) = ROI1(3)+24;
            ROI1(4) = ROI1(4)-24;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 =ROI1;
        elseif trialtype(rep_ind) == 4 %moved
            ROI1 = ROIs{img}{1};%original location
            ROI1(1) = ROI1(1)-24;
            ROI1(2) = ROI1(2)+24;
            ROI1(3) = ROI1(3)+24;
            ROI1(4) = ROI1(4)-24;
            ROI1(ROI1 < 1) = 1;
            ROI1(ROI1 > imageX) = imageX;
            if ROI1(3) > imageY; ROI1(3) = imageY;end
            if ROI1(4) > imageY; ROI1(4) = imageY;end
            ROI2 = ROIs{img}{2};%new location after move
            ROI2(1) = ROI2(1)-24;
            ROI2(2) = ROI2(2)+24;
            ROI2(3) = ROI2(3)+24;
            ROI2(4) = ROI2(4)-24;
            ROI2(ROI2 < 1) = 1;
            ROI2(ROI2 > imageX) = imageX;
            if ROI2(3) > imageY; ROI2(3) = imageY;end
            if ROI2(4) > imageY; ROI2(4) = imageY;end
        end
        
%           %added commented out code below to ploty eye traces on
%         %images.
%         novx = 24*fixationstats{nov_ind}.XY(1,1:600)+400;
%         novy = 24*fixationstats{nov_ind}.XY(2,1:600)+300;
%         
%         repx = 24*fixationstats{rep_ind}.XY(1,1:600)+400;
%         repy = 24*fixationstats{rep_ind}.XY(2,1:600)+300;
%         
%         scm_image_dir = ['R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\'...
%             'Set ' num2str(sets(SET)) '\Set' num2str(sets(SET)) '\'];
%         rep_trialType = trialtype(rep_ind);
%         
%         figure
%         subplot(2,1,1)
%         imshow(imread([scm_image_dir num2str(which_image(nov_ind)) '.bmp']));
%         box off
%         hold on
%         plot(novx,imageY-novy);
%         plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
%             600-[ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
%         plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
%             600-[ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
%         hold off
%         
%         subplot(2,1,2)
%         if rep_trialType == 2 % familiar
%             imshow(imread([scm_image_dir num2str(which_image(rep_ind)) 'p.bmp']));
%         elseif rep_trialType ==3 %replaced
%             imshow(imread([scm_image_dir num2str(which_image(rep_ind)) 'r.bmp']));
%         elseif rep_trialType ==4 %moved
%             imshow(imread([scm_image_dir num2str(which_image(rep_ind)) 'm.bmp']));
%         end
%         box off
%         hold on
%         plot(repx,imageY-repy);
%         plot([ROI1(1) ROI1(2) ROI1(2) ROI1(1) ROI1(1)],...
%             600-[ROI1(3) ROI1(3) ROI1(4) ROI1(4) ROI1(3)],'g');
%         plot([ROI2(1) ROI2(2) ROI2(2) ROI2(1) ROI2(1)],...
%             600-[ROI2(3) ROI2(3) ROI2(4) ROI2(4) ROI2(3)],'r');
%         hold off
%         
        ROIfix2 = find(...
            repfixations(1,:) > ROI2(1) & repfixations(1,:) < ROI2(2) & ...
            repfixations(2,:) < ROI2(3) & repfixations(2,:) > ROI2(4));%ROIfix2 contains indices of repfixations corresponding to fixations that occur in ROI of the manipulated image
        ROIfix1 = find(...
            novfixations(1,:) > ROI1(1) & novfixations(1,:) < ROI1(2) & ...
            novfixations(2,:) < ROI1(3) & novfixations(2,:) > ROI1(4)); %ROIfix1 contains indices of novfixations corresponding to fixations that occur in ROI of novel image.
        if trialtype(rep_ind) ==4 %moved image. Data from the old ROI (where the object was before being moved)
            ROIfix2_oldROI = find(...
                (repfixations(1,:) > ROI1(1) & repfixations(1,:) < ROI1(2) & ...
                repfixations(2,:) < ROI1(3) & repfixations(2,:) > ROI1(4)));%ROIfix2_oldROI contains indices of repfixations corresponding to fixations that occur in the 'old' ROI of the 'moved' image
        else
            ROIfix2_oldROI = [];
        end
%         close
        if isempty(ROIfix1)%didnt look at original region, moving on to next image
            continue
        end
        
        %%%%Code for calculating avg time in ROI for each image type
        tempvec1=zeros(1,5000);%novel image
        for j=1:length(ROIfix1)
            tempvec1(novtimes(1,ROIfix1(j)):novtimes(2,ROIfix1(j)))=1;
        end
        
        tempvec2=zeros(1,5000);%manipulated image new ROI
        for j=1:length(ROIfix2)
            tempvec2(reptimes(1,ROIfix2(j)):reptimes(2,ROIfix2(j)))=1;
        end
        
        tempvec3=zeros(1,5000);%manipulated imaged old moved ROI
        for j=1:length(ROIfix2_oldROI)
            tempvec3(reptimes(1,ROIfix2_oldROI(j)):reptimes(2,ROIfix2_oldROI(j)))=1;
        end
        
        area_ROI1 = abs(ROI1(1)-ROI1(2))*abs(ROI1(3)-ROI1(4)); %area of ROI for novel image
        area_ROI2 = abs(ROI2(1)-ROI2(2))*abs(ROI2(3)-ROI2(4)); %for replaced, moved or familiar
        
        %for novel presentation
        numFixationsInROI{SET,1} = [numFixationsInROI{SET,1} length(ROIfix1)];
        area{SET,1} = [area{SET,1} area_ROI1];
        Time_ROI_timeWindow{SET,1} = sum(tempvec1);
        
        %for second presentation
        numFixationsInROI{SET,trialtype(rep_ind)} = ...
            [numFixationsInROI{SET,trialtype(rep_ind)} length(ROIfix2)+length(ROIfix2_oldROI)];
        if trialtype(rep_ind) == 4
            area{SET,trialtype(rep_ind)} = [area{SET,trialtype(rep_ind)} area_ROI2+area_ROI1];
        else
            area{SET,trialtype(rep_ind)} = [area{SET,trialtype(rep_ind)} area_ROI2];
        end
        Time_ROI_timeWindow{SET,trialtype(rep_ind)} =  [...
            Time_ROI_timeWindow{SET,trialtype(rep_ind)} sum(tempvec2+tempvec3)];
        
        prop_fixROI_timeWindow = cell(length(cortex_files),4);
        
        
        all_time(1,:) = all_time(1,:) +tempvec1;
        all_time(trialtype(rep_ind),:) =  all_time(trialtype(rep_ind),:)+tempvec2+tempvec3;
        all_time(5,:) = all_time(5,:) + tempvec3;%for new moved only
    end
    
    clrs = 'brmgg';
    line = {'-','-','-','-','--'};
    image_count = cellfun(@length,numFixationsInROI(SET,:));%num images in each set
    image_count = [image_count image_count(end)];
    figure
    hold on
    for r = 1:size(all_time,1)
        plot(100*filtfilt(1/200*ones(1,200),1,all_time(r,:)/image_count(r)),...
            [clrs(r) line{r}]);
    end
    hold off
    xlabel('Time (ms)')
    ylabel('% of time in ROI')
    legend('Novel','Repeat','Replaced','Moved (new+old)','Moved (new)');
    title(['Set ' num2str(sets(SET)) ' ' num2str(imgdur(SET)) ' sec presentatoin  with ' ...
        num2str(img_seperation(SET)) 'images between image 1st and 2nd presentation'])
end
%%
SET =2 ; 
figure
bar(cellfun(@mean,numFixationsInROI(SET,:))/10)
set(gca,'XtickLabel',{'Novel','Repeat','Replaced','Moved'})
ylabel('Proportion of 1st 10 Fixations in ROI')
   title(['Set ' num2str(sets(SET)) ' ' num2str(imgdur(SET)) ' sec presentatoin  with ' ...
        num2str(img_seperation(SET)) 'images between image 1st and 2nd presentation'])
    %%
      image_count = cellfun(@length,numFixationsInROI(SET,:));%num images in each set
    image_count = [image_count image_count(end)];
    average_time = [];
    for r = 1:size(all_time,1)
        average_time(r) = 100*mean(all_time(r,1:3000)/image_count(r));
    end