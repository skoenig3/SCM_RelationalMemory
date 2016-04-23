% Code written to Analyze Human SCM behavioral Performance
% Written by Seth Konig 3/27/15
img_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\';
b2b_dir = 'B2B_responses\';%responses to back-2-back SCM
list12_dir = 'List12_responses\';

%---Analyze Back to Back Responses---%
a = what([img_dir b2b_dir]);
matfiles = a.mat;

b2b_numtype = NaN(size(matfiles,1),4);%number of response for each type
b2b_type_correct = NaN(size(matfiles,1),3);%percent correct by type
b2b_numNS = NaN(size(matfiles,1),4);%percent of response that were not sure for each type

b2b_by_img = zeros(15,36); %number of correct response by image
b2b_type = zeros(15,36);%image type
b2b_num_response = zeros(15,36);%total number of response per image

for file = 1:size(matfiles,1)
    fname = matfiles{file};
    load([img_dir b2b_dir fname])
    underscore = strfind(fname,'_');
    setname = fname(underscore(1)+1:underscore(2)-1);%get the image set
    folder_name = [img_dir setname]; %image folder
    
    %---get the name of images in the foldler---%
    image_list = ls([folder_name '\','*.bmp']);
    order = NaN(2,36);
    for il = 1:size(image_list,1);
        num  = sscanf(image_list(il,:), '%d');
        if num < 10
            str  = sprintf(image_list(il,2:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        else
            str  = sprintf(image_list(il,3:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        end
        order(row,num) = il;
    end
    
    %if organizing back-to-back
    image_names = [];
    for img = 1:36;
        %strtrim removes extra spaces if they exist in the string names
        image_names = [image_names {strtrim(image_list(order(1,img),:))} ...
            {strtrim(image_list(order(2,img),:))}];
    end
    
    image_type = NaN(1,72);%1 novel, 2 repeat, 3 replaced, 4 moved
    image_type(1:2:end) = 1;%novel images
    for im =2:2:72
        if ~isempty(strfind(image_names{im}(1:end-4),'p'))%repeat
            image_type(im) = 2;
        elseif ~isempty(strfind(image_names{im}(1:end-4),'r'))%replaced
            image_type(im) = 3;
        elseif ~isempty(strfind(image_names{im}(1:end-4),'m'))%moved
            image_type(im) = 4;
        else
            error('Image type unknown')
        end
    end
    
    %---Compare users responses to actual image types--%
    image_type(1:2:end) = [];%remove novels
    Responses(1:2:end) =[];%removed novels for responses, they were never probed
    
    
    %calculate the number of response for each type
    for type = 2:5
        b2b_numtype(file,type-1) = sum(Responses == type); %number of response for each type
    end
    
    
    %store data across sets and people
    setnum = str2double(setname(4:end));
    correct_trials = find(Responses == image_type);
    b2b_by_img(setnum,correct_trials) = b2b_by_img(setnum,correct_trials)+1;
    b2b_type(setnum,:) = image_type;%image type, going to write over
    b2b_num_response(setnum,:) = b2b_num_response(setnum,:)+1;%total number of response per image
    
    %calculate the number of correct response for each type
    for type = 2:4
        this_type = find(image_type == type);
        these_responses = Responses(this_type);
        b2b_type_correct(file,type-1) = 100*sum(these_responses == type)/12; %percent correct by type
        b2b_numNS(file,type-1) = 100*sum(these_responses == 5)/12;%percent of response that were not sure for each type
    end
end

figure
hold on
bar(mean(b2b_numtype))
errorb(mean(b2b_numtype),std(b2b_numtype)./sqrt(size(b2b_numtype,1)))
set(gca,'Xtick',[1:4])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved','Not Sure'})
ylabel('Average Number of Response')
title('B2B Number of Response by Response Type')

figure
hold on
bar(mean(b2b_type_correct))
errorb(mean(b2b_type_correct),std(b2b_type_correct)./sqrt(size(b2b_type_correct,1)))
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Correct')
title('B2B Percentage of Correct Response by Type')

figure
hold on
bar(mean(b2b_numNS(:,1:3)))
errorb(mean(b2b_numNS(:,1:3)),std(b2b_numNS(:,1:3))./sqrt(size(b2b_numNS(:,1:3),1)))
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Responses')
title('B2B Percentage of Response that were Not Sure by Type')

%---Analyze List 12 style Responses---%
a = what([img_dir list12_dir]);
matfiles = a.mat;

List12_numtype = NaN(size(matfiles,1),4);%number of response for each type
List12_type_correct = NaN(size(matfiles,1),3);%percent correct by type
List12_numNS = NaN(size(matfiles,1),4);%percent of response that were not sure for each type

list12_by_img = zeros(15,36); %number of correct response by image
list12_type = zeros(15,36);%image type
list12_num_response = zeros(15,36);%total number of response per image


for file = 1:size(matfiles,1)
    fname = matfiles{file};
    load([img_dir list12_dir fname])
    underscore = strfind(fname,'_');
    setname = fname(underscore(1)+1:underscore(2)-1);%get the image set
    folder_name = [img_dir setname]; %image folder
    
    %---get the name of images in the foldler---%
    image_list = ls([folder_name '\','*.bmp']);
    order = NaN(2,36);
    for il = 1:size(image_list,1);
        num  = sscanf(image_list(il,:), '%d');
        if num < 10
            str  = sprintf(image_list(il,2:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        else
            str  = sprintf(image_list(il,3:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        end
        order(row,num) = il;
    end
    
    %---get the name of images in the foldler---%
    image_list = ls([folder_name '\','*.bmp']);
    order = NaN(2,36);
    for il = 1:size(image_list,1);
        num  = sscanf(image_list(il,:), '%d');
        if num < 10
            str  = sprintf(image_list(il,2:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        else
            str  = sprintf(image_list(il,3:end), '%s');
            if strcmpi(str(1),'.')
                row = 1;%novel image
            else
                row = 2;%manipulated or repeat image
            end
        end
        order(row,num) = il;
    end
    
    %if organizing List 12 style
    image_type = []; %where we're going to store the data, row1 type, row 2 guess type
    image_names = [];
    for block = 1:3
        for block_img = 1:12
            img = 12*(block-1)+block_img;
            image_names = [image_names {strtrim(image_list(order(1,img),:))}];
        end
        for block_img = 1:12;
            img = 12*(block-1)+block_img;
            image_names = [image_names {strtrim(image_list(order(2,img),:))}];
        end
    end
    
    image_type = NaN(1,72);%1 novel, 2 repeat, 3 replaced, 4 moved
    image_type(1:2:end) = 1;%novel images
    for im =1:72
        if ~isempty(strfind(image_names{im}(1:end-4),'p'))%repeat
            image_type(im) = 2;
        elseif ~isempty(strfind(image_names{im}(1:end-4),'r'))%replaced
            image_type(im) = 3;
        elseif ~isempty(strfind(image_names{im}(1:end-4),'m'))%moved
            image_type(im) = 4;
        else
            image_type(im) = 1;%novel
        end
    end
    
    %calculate the number of response for each type
    for type = 2:5
        List12_numtype(file,type-1) = sum(Responses == type); %number of response for each type
    end
    
    
    %store data across sets and people
    image_type2 = image_type(image_type ~= 1);
    Responses2 = Responses(image_type ~= 1);
    setnum = str2double(setname(4:end));
    correct_trials = find(Responses2 == image_type2);
    list12_by_img(setnum,correct_trials) = list12_by_img(setnum,correct_trials)+1;
    list12_type(setnum,:) = image_type2;%image type, going to write over
    list12_num_response(setnum,:) = list12_num_response(setnum,:)+1;%total number of response per image
    
    %calculate the number of correct response for each type
    for type = 2:4
        this_type = find(image_type == type);
        these_responses = Responses(this_type);
        List12_type_correct(file,type-1) = 100*sum(these_responses == type)/12; %percent correct by type
        List12_numNS(file,type-1) = 100*sum(these_responses == 5)/12;%percent of response that were not sure for each type
    end
end

figure
hold on
bar(mean(List12_numtype))
errorb(mean(List12_numtype),std(List12_numtype)./sqrt(size(List12_numtype,1)))
set(gca,'Xtick',[1:4])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved','Not Sure'})
ylabel('Average Number of Response')
title('List12 Number of Response by Response Type')

figure
hold on
bar(mean(List12_type_correct))
errorb(mean(List12_type_correct),std(List12_type_correct)./sqrt(size(List12_type_correct,1)))
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Correct')
title('List12 Percentage of Correct Response by Type')

figure
hold on
bar(mean(List12_numNS(:,1:3)))
errorb(mean(List12_numNS(:,1:3)),std(List12_numNS(:,1:3))./sqrt(size(List12_numNS(:,1:3),1)))
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Responses')
title('List12 Percentage of Response that were Not Sure by Type')

%%
%---Plot both List12 and B2B results on same plot---%
numtype_means = [mean(List12_numtype)' mean(b2b_numtype)'];
numtype_sems =  [std(List12_numtype)'./sqrt(size(List12_numtype,1)) ...
    std(b2b_numtype)'./sqrt(size(b2b_numtype,1))];
[~,numtype_pvals(1)] = ttest2(List12_numtype(:,1),b2b_numtype(:,1));
[~,numtype_pvals(2)] = ttest2(List12_numtype(:,2),b2b_numtype(:,2));
[~,numtype_pvals(3)] = ttest2(List12_numtype(:,3),b2b_numtype(:,3));
[~,numtype_pvals(4)] = ttest2(List12_numtype(:,4),b2b_numtype(:,4));

figure
hold on
errorb(numtype_means,numtype_sems)
plot([1 3],[19 10],'*k')
hold off
set(gca,'Xtick',[1:4])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved','Not Sure'})
ylabel('Average Number of Response')
title('Number of Response by Response Type')
legend('List12','B2B','Location','NorthEastOutside')

type_correct_means = [mean(List12_type_correct)' mean(b2b_type_correct)'];
type_correct_sem = [std(List12_type_correct)'./sqrt(size(List12_type_correct,1))...
    std(b2b_type_correct)'./sqrt(size(b2b_type_correct,1))];
[~,type_correct_pvals(1)] = ttest2(List12_type_correct(:,1),b2b_type_correct(:,1));
[~,type_correct_pvals(2)] = ttest2(List12_type_correct(:,2),b2b_type_correct(:,2));
[~,type_correct_pvals(3)] = ttest2(List12_type_correct(:,3),b2b_type_correct(:,3));

figure
hold on
errorb(type_correct_means,type_correct_sem)
plot([1 3],[95 70],'*k')
hold off
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Correct')
title('Percentage of Correct Response by Type')
legend('List12','B2B','Location','NorthEastOutside')

numNs_means = [mean(List12_numNS(:,1:3))' mean(b2b_numNS(:,1:3))'];
numNs_sems = [std(List12_numNS(:,1:3))'./sqrt(size(List12_numNS(:,1:3),1))...
    std(b2b_numNS(:,1:3))'./sqrt(size(b2b_numNS(:,1:3),1))];
[~,Ns_pvals(1)] = ttest2(List12_numNS(:,1),b2b_numNS(:,1));
[~,Ns_pvals(2)] = ttest2(List12_numNS(:,2),b2b_numNS(:,2));
[~,Ns_pvals(3)] = ttest2(List12_numNS(:,3),b2b_numNS(:,3));

figure
errorb(numNs_means,numNs_sems)
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Responses')
legend('List12','B2B')
legend('List12','B2B','Location','NorthEastOutside')
title('Percentage of Response that were Not Sure by Type')

%%
%determine if the images everyone got correct on b2b were harder during
%list12 style
b2b_percent_correct = b2b_by_img./b2b_num_response;
list12_percent_correct = list12_by_img./list12_num_response;

means = NaN(1,3);
ns = NaN(1,3);
for type = 2:4
    img_ind = find(b2b_percent_correct == 1 & b2b_type == type);
    means(type-1) = mean(list12_percent_correct(img_ind)); 
    ns(type-1) = length(img_ind);
end

figure
bar(100*means)
set(gca,'Xtick',[1:3])
set(gca,'XtickLabel',{'Repeat','Replaced','Moved'})
ylabel('% Correct')
title('Probability of Correct Response on List12 style given all subjects got type correct on B2B')
%% Write data to excel spread sheet
type_letter= 'bprm';
array = cell(1,15);
for set = 1:15
    sub_array = cell(38,2);
    sub_array{1,1} = ['Set ' num2str(set)];
    sub_array{2,1} = 'Image #';
    sub_array{2,2} = '% Correct'; 
    for img = 1:36
        sub_array{img+2,1} = [num2str(img) type_letter(b2b_type(set,img))];
        sub_array{img+2,2} = 100*b2b_percent_correct(set,img);
    end
    array{set} = sub_array; 
end

filename ='Human_SCM_B2B_Results';
for set = 1:15
    xlswrite(filename,array{set},set)
end