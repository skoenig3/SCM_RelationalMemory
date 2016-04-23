%written by Seth Konig ok 8/13/15
%code takes a bunch of repeated, moved, and replaced images and sorts them
%into image folders. Also makes appropriate item files. If images pools
%have not changed then code can exactly replicate sets if images are lost

main_img_dir = 'C:\Users\seth.koenig\Desktop\';
set_dir = 'C:\Users\seth.koenig\Documents\MATLAB\';

%image pools by type
repeat_dir = [main_img_dir 'All repeats\'];
replaced_dir = [main_img_dir 'All replaced\'];
moved_dir = [main_img_dir 'All moved\'];

%get the image names in each folder
repeat_list = ls(repeat_dir);
replaced_list = ls(replaced_dir);
moved_list = ls(moved_dir);

%---Check to make sure that all types of images have a novel and their respective repeat/manipulation version---%
%also get the novel list too
%assumes images are not in the wrong folder
novel_repeat_list_index = NaN(1,480);
imgnum = 1;
for l = 1:length(repeat_list)
    bmp = strfind(repeat_list(l,:),'.bmp');
    if ~isempty(bmp)
        if ~strcmp(repeat_list(l,bmp-1),'p') %then novel image
            if ~exist([repeat_dir repeat_list(l,1:bmp-1) 'p.bmp'],'file')
                disp(['Repeat pairing missing for ' repeat_list(l,:)])
            else
                novel_repeat_list_index(imgnum) = l;
                imgnum = imgnum+1;
            end
        end
    end
end

novel_moved_list_index = NaN(1,480);
imgnum = 1;
for l = 1:length(moved_list)
    bmp = strfind(moved_list(l,:),'.bmp');
    if ~isempty(bmp)
        if ~strcmp(moved_list(l,bmp-1),'m') %then novel image
            if ~exist([moved_dir moved_list(l,1:bmp-1) 'm.bmp'],'file')
                disp(['Moved pairing missing for ' moved_list(l,:)])
            else
                novel_moved_list_index(imgnum) = l;
                imgnum = imgnum+1;
            end
        end
    end
end

novel_replaced_list_index = NaN(1,480);
imgnum = 1;
for l = 1:length(replaced_list)
    bmp = strfind(replaced_list(l,:),'.bmp');
    if ~isempty(bmp)
        if ~strcmp(replaced_list(l,bmp-1),'r') %then novel image
            if ~exist([replaced_dir replaced_list(l,1:bmp-1) 'r.bmp'],'file')
                disp(['Replaced pairing missing for ' replaced_list(l,:)])
            else
                novel_replaced_list_index(imgnum) = l;
                imgnum = imgnum+1;
            end
        end
    end
end

if length(novel_repeat_list_index) < 480
    error('Not enough repeat images to create 40 sets')
end
if length(novel_replaced_list_index) < 480
    error('Not enough repeat images to create 40 sets')
end
if length(novel_moved_list_index) < 480
    error('Not enough repeat images to create 40 sets')
end
%%
%---Sort Images into Sets---%
rand('seed',81315) %seed so can redo everything the same way

%randomize the image orders
randind = randperm(length(novel_repeat_list_index));
novel_repeat_list_index = novel_repeat_list_index(randind);

randind = randperm(length(novel_replaced_list_index));
novel_replaced_list_index = novel_replaced_list_index(randind);

randind = randperm(length(novel_moved_list_index));
novel_moved_list_index = novel_moved_list_index(randind);

repeat_img_ind = 1;
replaced_img_ind = 1;
moved_img_ind = 1;
for set = 1:40
    if set < 10
        new_set_dir = [set_dir 'SCMRM0' num2str(set) '\'];
    else
        new_set_dir = [set_dir 'SCMRM' num2str(set) '\'];
    end
    mkdir(new_set_dir)
    
    for block = 1:3
        repeat_order = [1 1 2 2];%1/2 of images should have presented manipulated image presented before original
        replaced_order = [1 1 2 2];%1/2 of images should have presented manipulated image presented before original
        moved_order = [1 1 2 2];%1/2 of images should have presented manipulated image presented before original
        repeat_order = repeat_order(randperm(4));
        replaced_order = replaced_order(randperm(4));
        moved_order = moved_order(randperm(4));
        
        type = [2 2 2 2 3 3 3 3 4 4 4 4]; %2 for repeats, 3 for replacd, 4 for moved
        type = type(randperm(length(type))); %rnadomize order for each block
        repeat_ind = find(type == 2);
        replaced_ind = find(type == 3);
        moved_ind = find(type == 4);
        base_image_num = 12*(block-1);
        
        for im = 1:4
            if set < 10
                if repeat_ind(im)+base_image_num < 10
                    new_img_name = ['S0' num2str(set) 'I0' num2str(repeat_ind(im)+base_image_num)];
                else
                    new_img_name = ['S0' num2str(set) 'I' num2str(repeat_ind(im)+base_image_num)];
                end
            else
                if repeat_ind(im)+base_image_num < 10
                    new_img_name = ['S' num2str(set) 'I0' num2str(repeat_ind(im)+base_image_num)];
                else
                    new_img_name = ['S' num2str(set) 'I' num2str(repeat_ind(im)+base_image_num)];
                end
            end
            bmp = strfind(repeat_list(novel_repeat_list_index(repeat_img_ind),:),'.bmp');
            old_img_name = repeat_list(novel_repeat_list_index(repeat_img_ind),1:bmp-1);
            if repeat_order(im) == 1  %original first followed by manipulated
                copyfile([repeat_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([repeat_dir old_img_name 'p.bmp'],...
                    [new_set_dir new_img_name 'p.bmp']);
            else  %swap so manipulated first followed by original
                copyfile([repeat_dir old_img_name 'p.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([repeat_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name 'p.bmp']);
            end
            repeat_img_ind = repeat_img_ind+1;
        end
        
        for im = 1:4
            if set < 10
                if replaced_ind(im)+base_image_num < 10
                    new_img_name = ['S0' num2str(set) 'I0' num2str(replaced_ind(im)+base_image_num)];
                else
                    new_img_name = ['S0' num2str(set) 'I' num2str(replaced_ind(im)+base_image_num)];
                end
            else
                if replaced_ind(im)+base_image_num < 10
                    new_img_name = ['S' num2str(set) 'I0' num2str(replaced_ind(im)+base_image_num)];
                else
                    new_img_name = ['S' num2str(set) 'I' num2str(replaced_ind(im)+base_image_num)];
                end
            end
            bmp = strfind(replaced_list(novel_replaced_list_index(replaced_img_ind),:),'.bmp');
            old_img_name = replaced_list(novel_replaced_list_index(replaced_img_ind),1:bmp-1);
            if replaced_order(im) == 1 %original first followed by manipulated
                copyfile([replaced_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([replaced_dir old_img_name 'r.bmp'],...
                    [new_set_dir new_img_name 'r.bmp']);
            else %swap so manipulated first followed by original
                copyfile([replaced_dir old_img_name 'r.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([replaced_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name 'r.bmp']);
            end
            replaced_img_ind = replaced_img_ind+1;
        end
        
        for im = 1:4
            if set < 10
                if moved_ind(im)+base_image_num < 10
                    new_img_name = ['S0' num2str(set) 'I0' num2str(moved_ind(im)+base_image_num)];
                else
                    new_img_name = ['S0' num2str(set) 'I' num2str(moved_ind(im)+base_image_num)];
                end
            else
                if moved_ind(im)+base_image_num < 10
                    new_img_name = ['S' num2str(set) 'I0' num2str(moved_ind(im)+base_image_num)];
                else
                    new_img_name = ['S' num2str(set) 'I' num2str(moved_ind(im)+base_image_num)];
                end
            end
            bmp = strfind(moved_list(novel_moved_list_index(moved_img_ind),:),'.bmp');
            old_img_name = moved_list(novel_moved_list_index(moved_img_ind),1:bmp-1);
            if moved_order(im) == 1  %original first followed by manipulated
                copyfile([moved_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([moved_dir old_img_name 'm.bmp'],...
                    [new_set_dir new_img_name 'm.bmp']);
            else  %swap so manipulated first followed by original
                copyfile([moved_dir old_img_name 'm.bmp'],...
                    [new_set_dir new_img_name '.bmp']);
                copyfile([moved_dir old_img_name '.bmp'],...
                    [new_set_dir new_img_name 'm.bmp']);
            end
            moved_img_ind = moved_img_ind+1;
        end
    end
end

%% Make item files

main_img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\';

spacingx = {'-12.0','-6.00',' 0.00',' 6.00',' 12.0'};
spacingy = {'-8.00','-4.00',' 0.00',' 4.00',' 8.00'};
ngrid = length(spacingx);

numimages = 72; %36 novel as well as 12 repeats, 12 moved, and 12 replaced
for setnum = 1:40
    if setnum < 10
        set = ['SCMRM0' num2str(setnum)];
    else
        set = ['SCMRM' num2str(setnum)];
    end
    
    fid = fopen([set '.itm'],'w+');
    
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
    line2 =' -4   14      1    0.00    0.00      0   0.50  0.50  0.00              255 255 255 x';
    line3 =' -3   14      1    0.00    0.00      0   0.50  0.50  0.00               75  75  75 x';
    line4 =' -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         10  10  10 x';
    line5 =' -1    1      0    0.00    0.00      0   1.00  1.00  0.00              200 200 200 x';
    line6 ='  0    1      1    0.00    0.00      0   0.15  0.15  0.00              150 150 150 x';
    
    
    for line = 1:6
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    for i = 1:length(spacingx)
        for ii = 1:length(spacingy)
            num1 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii);
            num3 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii+1);
            num2 = num2str(0+3*ngrid*(i-1)+2*(ii-1)+ii+2);
            if str2double(num1) < 10
                header1 = ['  ' num1];
            else
                header1 = [' ' num1];
            end
            if str2double(num2) < 10
                header2 = ['  ' num2];
            else
                header2 = [' ' num2];
            end
            if str2double(num3) < 10
                header3 = ['  ' num3];
            else
                header3 = [' ' num3];
            end
            str1 = [header1 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        150 150 150 x' '\r\n'];
            str3 = [header3 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        150 150 150 x' '\r\n'];
            str2 = [header2 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50        175 175 130 x' '\r\n'];
            fprintf(fid,str1);
            fprintf(fid,str3);
            fprintf(fid,str2);
        end
    end
    
    if setnum < 10
        setzero = '0';
    else
        setzero = '';
    end
    
    
    if setnum < 10
        
        image_list = ls([main_img_dir 'SCM0' num2str(setnum) '\','*.bmp']);
    else
        image_list = ls([main_img_dir 'SCM' num2str(setnum) '\','*.bmp']);
    end
    
    order = NaN(2,36);
    for il = 1:size(image_list,1);
        num  = sscanf(image_list(il,5:6), '%d');
        str  = sprintf(image_list(il,7:end), '%s');
        
        if strcmpi(str(1),'.')
            row = 1;%novel image
        else
            row = 2;%manipulated or repeat image
        end
        order(row,num) = il;
    end
    
    image_names = cell(1,72);
    for img = 1:36;
        %strtrim removes extra spaces if they exist in the string names
        image_names{2*img-1} = strtrim(image_list(order(1,img),:)); %1st presentation
        image_names{2*img} = strtrim(image_list(order(2,img),:));%2nd presentation
    end
    
    
    for i = 1:numimages;
        if i < 10
            imgzero = '0';
        else
            imgzero = '';
        end
        if 75+i < 100
            str = [' ' num2str(75+i) '    8           0.00    0.00      0                                  75  75  75 x   C:\\'...
                'SCM' setzero num2str(setnum)  '\\' image_names{i} '\r\n'];
        else
            str = ['' num2str(75+i) '    8           0.00    0.00      0                                  75  75  75 x   C:\\' ...
                'SCM' setzero num2str(setnum)  '\\' image_names{i} '\r\n'];
        end
        fprintf(fid,str);
    end
    fclose(fid);
end