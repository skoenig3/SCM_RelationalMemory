function get_SCM_RM_eyedat(cortexfile,setnum)
%written 3/1/15 by Seth Konig
%adopted from get_List12M_eyedat
%Imports cortex data for Scene Manipulation task for the Relational Memory project,
% and then runs Cluster Fix to %detect fixations and saccades.
%
% Inputs:
%   1) cortexfile: cortex file containing behavioral data for task
%   2) setnum: set number assocatied with cortex file e.g. 7 for
%   List12S7.itm
%
% Output:
%   1)  Matlab file with the save fixations and saccade times as well the
%   behavioral data for that session

samprate = 5;%number of ms between samples for ISCAN i.e. 200 Hz

if strcmpi(cortexfile(1:2),'WR')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Wilbur\' cortexfile];
elseif strcmpi(cortexfile(1:2),'PW')
    cortexfile = ['R:\Buffalo Lab\Cortex Data\Vivian\' cortexfile];
end

% for Wilbur Pilot analysis
% ITMFile = ['R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\SCM Picture Sets - SR\'...
%     'List12S' num2str(setnum) '.itm'];
% CNDFile = ['R:\Buffalo Lab\eblab\Cortex Programs\Scene Manipulation\List 36 5 sec\' ...
%     'List12M' num2str(setnum) '.cnd']; %same cnd file as List12M task

if setnum < 10
    ITMFile = ['R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\'...
        'Item Files\SCMRM0' num2str(setnum) '.itm'];
else
    ITMFile = ['R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\'...
        'Item Files\SCMRM0' num2str(setnum) '.itm'];
end
CNDFile = 'R:\Buffalo Lab\eblab\Cortex Programs\SCM Relational Memory Project\scmRM.cnd';
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

itmlist = zeros(size(cndfil,1)-1,1); %which item is associated with which condition
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

%for color change files we want to know the position of the item in dva
spacex = [-12,-6,0,6,12];
spacey = [-8,-4,0,4,8];
ind_spacex = spacex;
ind_spacey = spacey;

which_image=NaN(1,72);%which item is assocaited with which image
%went cautious route  but should  be [1:12 1:12 13:24 13:24 25:36 25:36]
trialtype = NaN(1,72);%1 NOVEL, 2 FAMIILIAR, 3 REPLACED, 4 MOVED
for itm = 76:147
    str = textscan(itmfil(itm+6,:),'%s');
    imgstr = str{1}(end);
    slash = strfind(imgstr{1},'\');
    period = strfind(imgstr{1},'.');
    imgname = textscan(imgstr{1}(slash(end)+1:period-1),'%d');
    which_image(itm-75) = imgname{1};
    
    if  ~isempty(strfind(imgstr{1}(slash(end)+1:period-1),'p'))
        trialtype(itm-75) = 2; %familiar
    elseif ~isempty(strfind(imgstr{1}(slash(end)+1:period-1),'r'))
        trialtype(itm-75) = 3; %replaced
    elseif  ~isempty(strfind(imgstr{1}(slash(end)+1:period-1),'m'))
        trialtype(itm-75) = 4; %moved
    else
        trialtype(itm-75) = 1;%novel
    end
end

[time_arr,event_arr,eog_arr,~,~,~]  = get_ALLdata(cortexfile);

%---Import Clrchng Eye data so can calibrate image trial data---%
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) <= 75 %clrchng trials
        if size(find(event_arr(:,rptlop) == 3)) ~=0 %rewarded so sucessful clrchng trials
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            if length( perbegind) > 1
                perbegind = perbegind(2);
                perendind = perendind(2);
            end
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            trialtypeind = find(event_arr(:,rptlop) <= 4);
            if event_arr(blknumind,rptlop) == 1 %skip the 1st block it's for offset corrections
                continue
            end
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
            end
        end
    end
end

clear cnd
numrpt = size(per,2);
cnd = zeros(1,numrpt);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end

%---For Calibration with Eye tracking data with cp2tform---%
x = cell(length(spacex),length(spacey));
y = cell(length(spacex),length(spacey));
control = NaN(length(cnd),2);
clr = ['rgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmkrgbmk'];
figure
hold on
for k = 1:length(cnd)
    C = textscan(itmfil(itmlist(cnd(k)-1000)+5,:),'%d');
    control(k,:) = C{1}(4:5)';
    
    xi = find(C{1}(4) == ind_spacex);
    yi = find(C{1}(5) == ind_spacey);
    eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
    evenind = eyeind(logical(~rem(eyeind,2)));
    oddind =  eyeind(logical(rem(eyeind,2)));
    x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
    y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
    plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
end
if iscell(cortexfile(end-9:end))
    title(['Calibration transformation for ' cortexfile{1}(end-9:end)])
else
    title(['Calibration transformation for ' cortexfile(end-9:end)])
end

%Test for errors%
count = zeros(length(spacey),length(spacex));
for xi = 1:length(spacex);
    for yi = 1:length(spacey);
        count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
    end
end
if any(count < 5);
    disp('Calibration trial analysis incomplete or error')
    disp('Check number of calibration pionts or task not finished')
end

clear meanx meany
for k=1:numel(x)
    xss = x{k};
    low = mean(xss)-std(xss);
    high = mean(xss)+std(xss);
    xss(xss < low) = [];
    xss(xss > high) = [];
    meanx(k)=median(xss);
end
for k=1:numel(y)
    yss = y{k};
    low = mean(yss)-std(yss);
    high = mean(yss)+std(yss);
    yss(yss < low) = [];
    yss(yss > high) = [];
    meany(k)=median(y{k});
end

%format the control (known from item file) locations of items
controlx = [];
controly = [];
for i = 1:length(spacex);
    for ii = 1:length(spacey);
        controly = [controly spacey(i)];
        controlx = [controlx spacex(ii)];
    end
end

%calculate the calibration function
tform = cp2tform([controlx' controly'], [meanx' meany'],'affine'); %can use polynomial 3/4 as well
tform.forward_fcn = tform.inverse_fcn;

newx = [];
newy = [];
figure
hold on
for i = 1:length(controlx);
    plot(controlx(i),controly(i),'r+')
    [x,y] = tformfwd(tform,meanx(i),meany(i));
    plot(x,y,'*b')
    newx(i) = x;
    newy(i) = y;
end
title(['Calibration transformation for ' cortexfile(end-9:end)])
xlim([-17.5 17.5])
ylim([-12.5 12.5])


%---Import Image Trial data---%
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
new_eog_arr=[];
for rptlop = 1:numrpt
    if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) >= 76
        if size(find(event_arr(:,rptlop) == 200)) ~=0 %200 is code for succesfully finishing a trial
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            if length(perbegind) > 1
                perbegind = perbegind(2);
                perendind = perendind(2);
            end
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            trialtypeind = find(event_arr(:,rptlop) >= 1 & event_arr(:,rptlop) <= 4);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                per(valrptcnt).event = rptlop;
                per(valrptcnt).trialtype = event_arr(trialtypeind,rptlop);
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
            end
        end
    end
end

%---get eye data for only when fixation cross or picture is displayed---%
eyedat = cell(1,length(per));
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    
    %this isn't how I do this part any more but for consistency sake
    %I'm going to copy exactly what I did before, it's not wrong! SKD
    %2/24/15
    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
    
    if picend > length(horeog)*5
        picend =  length(horeog)*5;
    end
    
    eyedat{trlop}(1,:) = horeog(round(picstart/samprate):floor(picend/samprate));
    eyedat{trlop}(2,:) = vrteog(round(picstart/samprate):floor(picend/samprate));
end

%---Recalibrate and automatically scale eye data---%
for eye = 1:length(eyedat)
    x = eyedat{eye}(1,:);
    y = eyedat{eye}(2,:);
    [x,y] = tformfwd(tform,x,y);
    eyedat{eye} = [x;y];
end

%---Run Cluster Fix to Detect fixations and saccades---%
fixationstats = ClusterFixation_Final(eyedat,5/1000);

%---Save Data to File---%
save(['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\' ...
    cortexfile(end-9:end-2) '_' cortexfile(end) '-fixation'],'fixationstats','per',...
    'setnum','which_image','trialtype')