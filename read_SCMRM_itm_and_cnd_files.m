function [itmlist,clrchng_locations,first_img_item,imgs,image_names,trialtype] = read_SCMRM_itm_and_cnd_files(itemfile)
% writteen by Seth Konig May, 2015. Modified from read_ListSQ_itm_and_cnd_files.mat
% Function imports item file and and grabs condition file to determine which items
% are associated with which condition (itmlist) since conditions are randomly
% organized.

ITMFile = ['\\towerexablox.wanprc.org\Buffalo\eblab\Cortex Programs\SCM Relational Memory Project\Item Files\' itemfile];
CNDFile = '\\towerexablox.wanprc.org\Buffalo\eblab\Cortex Programs\SCM Relational Memory Project\SCMRM.cnd';

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

first_img_item = [];
for i = 6:size(itmfil,1); %first 5 are hearder, background, etc
    str = textscan(itmfil(i,:),'%s');
    if ~isempty(strfind(str{1}{end},'.bmp')) && isempty(first_img_item)
        first_img_item  = str2num((str{1}{1}));
        break
    end
end

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

%probably a more efficient way of doing this since there's only 25 color
%change locations, but it may be easier for processing later
% also get which images were displayed. Should be the same across all items
% and condition sets but what the hey.
imgs = zeros(2,36);
image_names = cell(2,36);
clrchng_locations = cell(1,length(itmlist));
for cnd = 1:length(itmlist)
    if itmlist(cnd) < first_img_item %then it is a clrchng trial
        str = textscan(itmfil(itmlist(cnd)+6,:),'%d');
        clrchng_locations{cnd} =  double(str{1}(4:5)); %location in dva
    else
        str = textscan(itmfil(itmlist(cnd)+6,:),'%s');
        imgnum = str2double(str{1}{end}(14:15));
        if imgs(1,imgnum) == 0; %novel image
            imgs(1,imgnum) = cnd;
            image_names{1,imgnum} = str{1}{end}(10:end);
        else
            imgs(2,imgnum) = cnd;
            image_names{2,imgnum} = str{1}{end}(10:end);
        end
    end
end

trialtype = ones(2,36); %1 for novel
for i = 1:size(imgs,2);
    rep_name = image_names{2,i}(7);
    if strcmpi(rep_name,'p')
        trialtype(2,i) = 2;
    elseif strcmpi(rep_name,'r')
        trialtype(2,i) = 3;
    elseif strcmpi(rep_name,'m')
        trialtype(2,i) = 4;
    else
        error('image type unkown')
    end
end
end