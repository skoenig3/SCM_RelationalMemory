%written by Seth Konig...modified slightly from Make_SCM_Item_Files.m 5/12/15
% %---Rearrange images so that manipulations get presented uniformly---%

Set = 'Set02';
set_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\Image Sets\';
new_set_dir = [set_dir Set 'N\'];
mkdir(new_set_dir);

image_list = ls([set_dir Set '\']); 

%want to keep track of all the images and to make sure they all exist
image_index = zeros(2,36); %zero no image, 1 novel, 2 repeated, 3 replaced, 4 moved

%find novel images
image_index = 1;
for l = 1:size(image_list);
    period = strfind(image_list(l,:),'.bmp');
    if ~isempty(period) %so it's an image
        if ~isnan(image_list(l,period-1))
            image_index(1,image_number) = l;%label this image as found and novel
            image_number = image_number + 1;
        end
    end
end

if any(image_index(1,:) == 0)
    error('Not enought novel images present')
end
%% identify the type of 2nd presentation 
for novel = 1:36
    image_name = image_list(image_index(1,novel));
     period = strfind(image_name,'.bmp');
    if exist([set_dir Set '\']
    
    
%%
for l = 1:size(repeat_list);
    period = strfind(repeat_list(l,:),'.bmp');
    if ~isempty(period) %so it's an image
        image_number = str2double(repeat_list(l,1:period-2));
        image_type(2,image_number) = 2;%label this image as found and novel
    end
end

for l = 1:size(replaced_list);
    period = strfind(replaced_list(l,:),'.bmp');
    if ~isempty(period) %so it's an image
        image_number = str2double(replaced_list(l,1:period-2));
        image_type(2,image_number) = 3;%label this image as found and novel
    end
end

for l = 1:size(moved_list);
    period = strfind(moved_list(l,:),'.bmp');
    if ~isempty(period) %so it's an image
        image_number = str2double(moved_list(l,1:period-2));
        image_type(2,image_number) = 4;%label this image as found and novel
    end
end

repeats = find(image_type(2,:) == 2);
replaced = find(image_type(2,:) == 3);
moved = find(image_type(2,:) == 4);

last_image_num = 0;
for block = 1:3
    block_repeats = repeats(4*(block-1)+1:4*block);
    block_replaced = replaced(4*(block-1)+1:4*block);
    block_moved = moved(4*(block-1)+1:4*block);
    type = [2 2 2 2 3 3 3 3 4 4 4 4];
    type = type(randperm(length(type)));
    repeat_ind = find(type == 2);
    replaced_ind = find(type == 3);
    moved_ind = find(type == 4);
    base_image_num = 12*(block-1);
    
    for rep = 1:4
        copyfile([set_dir 'Original Images\' num2str(block_repeats(rep)) '.bmp'],...
            [new_set_dir 'Original Images\' num2str(repeat_ind(rep)+base_image_num) '.bmp']);%copy original image
        copyfile([set_dir 'Original Images\' num2str(block_repeats(rep)) '.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\' num2str(repeat_ind(rep)+base_image_num) '.bmp']);%copy original image into new set folder
        copyfile([set_dir 'Repeat\' num2str(block_repeats(rep)) 'p.bmp'],...
            [new_set_dir 'Repeat\' num2str(repeat_ind(rep)+base_image_num) 'p.bmp']);%copy repeated image
        copyfile([set_dir 'Repeat\' num2str(block_repeats(rep)) 'p.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\' num2str(repeat_ind(rep)+base_image_num) 'p.bmp']);%copy repeated image into new set folder
        copyfile([set_dir 'Repeat\ROI\' num2str(block_repeats(rep)) 'pr.bmp'],...
            [new_set_dir 'Repeat\ROI\' num2str(repeat_ind(rep)+base_image_num) 'pr.bmp']);%copy repeat image "ROI"
    end
    for repl = 1:4
        copyfile([set_dir 'Original Images\' num2str(block_replaced(repl)) '.bmp'],...
            [new_set_dir 'Original Images\' num2str(replaced_ind(repl)+base_image_num) '.bmp']);%copy original image
        copyfile([set_dir 'Original Images\' num2str(block_replaced(repl)) '.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\' num2str(replaced_ind(repl)+base_image_num) '.bmp']);%copy original into new set folder
        copyfile([set_dir 'Replaced\' num2str(block_replaced(repl)) 'r.bmp'],...
            [new_set_dir 'Replaced\' num2str(replaced_ind(repl)+base_image_num) 'r.bmp']);%copy replaceded image
        copyfile([set_dir 'Replaced\' num2str(block_replaced(repl)) 'r.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\' num2str(replaced_ind(repl)+base_image_num) 'r.bmp']);%copy replaceded image into new set folder
    end
    for mov = 1:4
        copyfile([set_dir 'Original Images\' num2str(block_moved(mov)) '.bmp'],...
            [new_set_dir 'Original Images\' num2str(moved_ind(mov)+base_image_num) '.bmp']);%copy original image
        copyfile([set_dir 'Original Images\' num2str(block_moved(mov)) '.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\'  num2str(moved_ind(mov)+base_image_num) '.bmp']);%copy original image into new set folder
        copyfile([set_dir 'Moved\' num2str(block_moved(mov)) 'm.bmp'],...
            [new_set_dir 'Moved\' num2str(moved_ind(mov)+base_image_num) 'm.bmp']);%copy moved image
        copyfile([set_dir 'Moved\' num2str(block_moved(mov)) 'm.bmp'],...
            [new_set_dir 'Set' num2str(Set) '\' num2str(moved_ind(mov)+base_image_num) 'm.bmp']);%copy moved image into new set folder
    end
end
end