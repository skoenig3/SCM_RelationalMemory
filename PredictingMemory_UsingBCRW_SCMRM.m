% Code written By Seth Konig to determine if certain maniuplations are
% harder than others to spot. Written January 25, 2016. Adpated from 
% BCRW_Manipulation_Predictive_Difficuly_Analysis

%%
%---[1] Make Salience maps for manipulated imagse---%
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\';
% 
% type = ['rm']; %object replaced and object moved
% 
% for SET =1:40;
%     
%     if SET < 10
%         setzero = '0';
%     else
%         setzero = '';
%     end
%     
%     
%     cd([scm_image_dir  'SCM' setzero (num2str(SET)) '\']);
%     
% 
%     for img = 1:36
%         
%         if img < 10
%             imgzero = '0';
%         else
%             imgzero = '';
%         end
%         
%         for t = 1:2
%             try
%                 getSalienceMap(['S' setzero num2str(SET) 'I' imgzero num2str(img) type(t) '.bmp'])
%             end
%         end
%     end
% end
% emailme('Finished creating saliencemaps')
%%
%---[2] Extract Viewing Behavior---%
% scm_eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
% 
% PLOTOPTIONS = 'none';
% imageX = 800;
% imageY = 600;
% SAMPRATE = 1;
% cd(scm_eyedata_dir)
% matfiles = what;
% eyedatafiles = [];
% for i = 1:length(matfiles.mat);
%     str = strfind(matfiles.mat{i},'fixation');
%     if ~isempty(str)
%         eyedatafiles = [eyedatafiles i];
%     end
% end
% for eyefile = eyedatafiles;
%     getViewingBehaviorSCMRM(matfiles.mat{eyefile},SAMPRATE,imageX,imageY,PLOTOPTIONS)
% end

%%
%---[3] Combine Viewing Behavior by Monkey---%
% scm_eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
% tags = {'PW','TT','TO','RR'};
% medianfix = NaN(20,length(tags));
% mediansac = NaN(20,length(tags));
% mediannumfix = NaN(20,length(tags));
% mediannumsac = NaN(20,length(tags));
% 
% cd(scm_eyedata_dir)
% 
% ind = ones(1,length(tags));
% matfiles = what;
% statfiles = [];
% for i = 1:length(matfiles.mat);
%     str = strfind(matfiles.mat{i},'ViewingBehavior');
%     if ~isempty(str)
%         for ii = 1:length(tags);
%             strt = strfind(matfiles.mat{i},tags{ii});
%             if ~isempty(strt)
%                 load(matfiles.mat{i},'fixduration','sacamplitude','avgsacprofile','avgfixprofile');
%                 mediannumfix(ind(ii),ii) = median(sum(~isnan(fixduration')));
%                 mediannumsac(ind(ii),ii) = median(sum(~isnan(sacamplitude')));
%                 medianfix(ind(ii),ii) = size(avgfixprofile,2);
%                 mediansac(ind(ii),ii) = size(avgsacprofile,2);
%                 ind(ii) = ind(ii)+1;
%             end
%         end
%     end
% end
% 
% %get the median number of fixations and saccades across all monkeys and sets
% mediannumfix = round(nanmedian(mediannumfix)); 
% mediannumsac = round(nanmedian(mediannumsac));
% %get the median "duration" of time warped data
% medianfix = round(nanmedian(medianfix));
% mediansac = round(nanmedian(mediansac));
% 
% allview = cell(1,length(tags));
% for i = 1:length(tags)
%     allview{i}.densitymap = zeros(600,800);
%     allview{i}.allfixations = [];
%     allview{i}.allsaccades = [];
%     allview{i}.persistence =[];
%     allview{i}.anglebtwfix = [];
%     allview{i}.sacangle_2fix = [];
%     allview{i}.distanceprofile = [];
%     allview{i}.distbtwnfix = [];
%     allview{i}.fixduration = [];
%     allview{i}.sacangle = [];
%     allview{i}.sacdist = [];
%     allview{i}.sacamplitude = [];
%     allview{i}.sacduration = [];
%     allview{i}.timebtwfix = [];
% end
% 
% cd(scm_eyedata_dir)
% 
% ind = ones(1,length(tags));
% matfiles = what;
% for i = 1:length(matfiles.mat);
%     str = strfind(matfiles.mat{i},'ViewingBehavior');
%     if ~isempty(str)
%         for ii = 1:length(tags);
%             strt = strfind(matfiles.mat{i},tags{ii});
%             if ~isempty(strt)
%                 load(matfiles.mat{i});
%                 if size(allfixations,2) == medianfix(ii);
%                     timewarp = 1:size(allfixations,2);
%                 else
%                     timewarp = round(linspace(1,size(allfixations,2),medianfix(ii)));
%                 end
%                 distanceprofile.fix = distanceprofile.fix(:,timewarp);
%                 persistence.fix = persistence.fix(:,timewarp);
%                 persistence.fix = persistence.fix(:,6:end-5); %first and last 5 are extended past eye movement
%                 distanceprofile.fix = distanceprofile.fix(:,6:end-5);  %first and last 5 are extended past eye movement
%                 allview{ii}.allfixations = [allview{ii}.allfixations;...
%                     allfixations(:,timewarp,:)];
%                 if size(allsaccades,2) == mediansac(ii);
%                     timewarp = 1:size(allsaccades,2);
%                 else
%                     timewarp = round(linspace(1,size(allsaccades,2),mediansac(ii)));
%                 end
%                 allview{ii}.allsaccades = [allview{ii}.allsaccades;...
%                     allsaccades(:,timewarp,:)];
%                 distanceprofile.sac = distanceprofile.sac(:,timewarp);
%                 distanceprofile.sac = distanceprofile.sac(:,6:end-5);  %first and last 5 are extended past eye movement
%                 persistence.sac = persistence.sac(:,timewarp);
%                 persistence.sac = persistence.sac(:,6:end-5); %first and last 5 are extended past eye movement
%                 try
%                 allview{ii}.persistence = [ allview{ii}.persistence;
%                     [persistence.sac persistence.fix]];
%                    allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
%                     [distanceprofile.sac distanceprofile.fix]];
%                 catch %probably a rounding error doesn't matter much for these if lose 1 sample
%                    disp('Inconsitent lengths') 
%                    rfix = size(persistence.fix,1);
%                    rsac = size(persistence.sac,1);
%                    if rsac > rfix
%                        dif = rsac-rfix; 
%                        allview{ii}.persistence = [ allview{ii}.persistence;
%                            [persistence.sac(1:end-dif,:) persistence.fix]];
%                        allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
%                            [distanceprofile.sac(1:end-dif,:) distanceprofile.fix]];
%                    else
%                        dif = rfix-rsac;
%                        allview{ii}.persistence = [ allview{ii}.persistence;
%                            [persistence.sac persistence.fix(1:end-dif,:)]];
%                        allview{ii}.distanceprofile = [allview{ii}.distanceprofile;
%                            [distanceprofile.sac distanceprofile.fix(1:end-dif,:)]];
%                    end
%                 end
%                 
%                 allview{ii}.anglebtwfix = [allview{ii}.anglebtwfix;anglebtwfix];
%                 allview{ii}.sacangle_2fix = [allview{ii}.sacangle_2fix;...
%                     sacangle_2fix];
%                 allview{ii}.densitymap = allview{ii}.densitymap+densitymap;
%              
%                 allview{ii}.distbtwnfix = [allview{ii}.distbtwnfix;distbtwnfix];
%                 allview{ii}.fixduration = [allview{ii}.fixduration;...
%                     fixduration];
%                 allview{ii}.sacangle = [allview{ii}.sacangle;sacangle];
%                 allview{ii}.sacdist = [allview{ii}.sacdist;sacdist];
%                 allview{ii}.sacamplitude = [allview{ii}.sacamplitude;sacamplitude];
%                 allview{ii}.sacduration = [allview{ii}.sacduration;sacduration];
%                 allview{ii}.timebtwfix = [allview{ii}.timebtwfix;...
%                     timebtwfix];
%                 allview{ii}.mediansac = mediansac(ii)-10; %reduce since the 1st 5 and last 5 are data from the other eye movement
%                 allview{ii}.medianfix = medianfix(ii)-10; %reduce since the 1st 5 and last 5 are data from the other eye movement
%                 allview{ii}.mediannumfix = mediannumfix(ii);
%                 allview{ii}.mediannumsac = mediannumsac(ii);
%             end
%         end
%     end
% end
% 
% clearvars -except image_sets allview tags
% 
% SAMPRATE = 1;
% n = (-180:180)*pi/180;
% variables = {'Dist','vel','accel','rot'};
% f = fspecial('gaussian',[256,256],24);
% graphnum = gcf;
% if graphnum == 1
%     graphnum = 0;
% end
% for i = 1:length(tags)
%     [allprobanglebtwfix] = hist(allview{i}.anglebtwfix(~isnan(allview{i}.anglebtwfix)),360);
%     allprobanglebtwfix = [allprobanglebtwfix(36:-1:1) allprobanglebtwfix allprobanglebtwfix(end:-1:end-36)];
%     allprobanglebtwfix = filtfilt(1/6*ones(1,6),1,allprobanglebtwfix);
%     allprobanglebtwfix = allprobanglebtwfix(37:end-37);
%     allprobanglebtwfix = allprobanglebtwfix/sum(allprobanglebtwfix);
%     allprobanglebtwfix = [allprobanglebtwfix allprobanglebtwfix(1)];
%     
%     figure(graphnum+1)
%     subtitle('Distribution of angles between fixations')
%     hax(1,i) = subplot(2,2,i);
%     polar(n,allprobanglebtwfix)
%     title(tags{i})
%     ph=findall(gca,'type','text');
%     set(ph,'fontweight','bold');
%     
%     [allprobsacangle] = hist(allview{i}.sacangle(~isnan(allview{i}.sacangle)),360);
%     allprobsacangle = [allprobsacangle(36:-1:1) allprobsacangle allprobsacangle(end:-1:end-36)];
%     allprobsacangle = filtfilt(1/6*ones(1,6),1,allprobsacangle);
%     allprobsacangle = allprobsacangle(37:end-37);
%     allprobsacangle = allprobsacangle/sum(allprobsacangle);
%     allprobsacangle = [allprobsacangle allprobsacangle(1)];
%     
%     figure(graphnum+2)
%     subtitle('Distribution of angles leaving a fixation')
%     hax(2,i) = subplot(2,2,i);
%     polar(n,allprobsacangle)
%     title(tags{i})
%     ph=findall(gca,'type','text');
%     set(ph,'fontweight','bold');
%     
%     %---Stats by fixation---%
%     allStatsbyfixation{i}.fixatoinspertrial = reshape(sum(~isnan(allview{i}.fixduration),2),[],36);
%     allStatsbyfixation{i}.meanfixationduration = SAMPRATE*nanmean(allview{i}.fixduration);
%     allStatsbyfixation{i}.stdfixationduration = SAMPRATE*nanstd(allview{i}.fixduration);
%     allStatsbyfixation{i}.numfix = sum(~isnan(allview{i}.fixduration));
%     allStatsbyfixation{i}.meansacdistance = nanmean(allview{i}.sacdist);
%     allStatsbyfixation{i}.stdsacdistance = nanstd(allview{i}.sacdist);
%     allStatsbyfixation{i}.numsacs = sum(~isnan(allview{i}.sacdist));
%     
%     figure(graphnum+3)
%     subtitle('Distribution of Fixation Durations')
%     hax(3,i) = subplot(2,2,i);
%     fixduration = allview{i}.fixduration(~isnan(allview{i}.fixduration))*SAMPRATE;
%     fixduration(fixduration > 500) = [];
%     hist(fixduration,95)
%     xlabel('Time (ms)')
%     title(tags{i})
%     
%     figure(graphnum+4)
%     subtitle('Distribution of Distances between Fixations')
%     hax(4,i) = subplot(2,2,i);
%     hist(allview{i}.distbtwnfix(~isnan(allview{i}.distbtwnfix)),100)
%     xlabel('Distance (Pixels)')
%     title(tags{i})
%     
%     figure(graphnum+5)
%     subtitle('Fixation Rate')
%     hax(5,i) = subplot(2,2,i);
%     hist(1000./(allview{i}.timebtwfix(~isnan(allview{i}.timebtwfix))*SAMPRATE),100)
%     xlabel('Hz')
%     title(tags{i})
%     
%     figure(graphnum+6)
%     subtitle('Distribution of Saccade Durations (Arc Lenght)')
%     hax(6,i) = subplot(2,2,i);
%     sacduration = allview{i}.sacduration(~isnan(allview{i}.sacduration))*SAMPRATE;
%     sacduration(sacduration > 150) = [];
%     hist(sacduration,28)
%     xlim([0 100])
%     xlabel('Time (ms)')
%     title(tags{i})
%     
%     
%     figure(graphnum+7)
%     subtitle('Distribution of Saccade Distances')
%     hax(7,i) = subplot(2,2,i);
%     sacdist = allview{i}.sacdist(~isnan(allview{i}.sacdist));
%     sacdist(sacdist > 800) = [];
%     hist(sacdist,100)
%     xlabel('Distance (Pixels)')
%     title(tags{i})
%     
%     figure(graphnum+8)
%     subtitle('Correlation between saccade duration and distance')
%     hax(8,i) = subplot(2,2,i);
%     plot(allview{i}.sacdist,allview{i}.sacduration*SAMPRATE,...
%         '*','markersize',2)
%     title(tags{i})
%     xlabel('Distance (pixels)')
%     ylabel('Duration (ms)')
%     xlim([0 1000])
%     ylim([0 200])
%     
%     figure(graphnum+9)
%     subtitle('Probability Distribution of Fixations')
%     hax(9,i) = subplot(2,2,i);
%     densitymap = allview{i}.densitymap;
%     densitymap = imfilter(densitymap,f);
%     densitymap = densitymap./sum(sum(densitymap));
%     imagesc(densitymap)
%     title(tags{i})
%     axis off
%     
%     figure(graphnum+10)
%     subtitle('Fixation Statistics by Trial Number, Fixation Number, or Saccade Number')
%     hax(10,i) = subplot(2,2,i);
%     title(tags{i})
%     hold all
%     errorbar(mean(allStatsbyfixation{i}.fixatoinspertrial),...
%         std(allStatsbyfixation{i}.fixatoinspertrial)/...
%         size(allStatsbyfixation{i}.fixatoinspertrial,1));
%     errorbar(allStatsbyfixation{i}.meanfixationduration,...
%         allStatsbyfixation{i}.stdfixationduration./sqrt(allStatsbyfixation{i}.numfix))
%     xl = find(allStatsbyfixation{i}.numfix < 50);
%     xlim([0 xl(1)])
%     ylim([0 400])
%     errorbar(allStatsbyfixation{i}.meansacdistance,...
%         allStatsbyfixation{i}.stdsacdistance./sqrt(allStatsbyfixation{i}.numsacs))
%     if i == 1;
%         ylabel('Number of Fixations, Fixation Duration (ms),Saccade Distance (Pixels)')
%         set(get(gca,'YLabel'),'Position',[-7 -0.3 0])
%     end
%     if i == 2
%         legend('# Fixations','Fixation Duration','Saccade Distance','Location',...
%             'NorthEastOutside')
%     end
%     if i > 2
%         xlabel('Trial Number, Fixation Number, or Saccade Number')
%     end
%     
%     avgfixation= mean(allview{i}.allfixations,1);
%     fixlen = size(avgfixation,2);
%     avgfixprofile = zeros(size(avgfixation));
%     for ii = 1:size(avgfixation,3);
%         avgfixprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,ii));
%         avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii) - min(avgfixprofile(:,:,ii));
%         avgfixprofile(:,:,ii) = avgfixprofile(:,:,ii)/max(avgfixprofile(:,:,ii));
%     end
%     avgsaccade= mean(allview{i}.allsaccades,1);
%     saclen = size(avgsaccade,2);
%     avgsacprofile = zeros(size(avgsaccade));
%     for ii = 1:size(avgsaccade,3);
%         avgsacprofile(:,:,ii) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,ii));
%         avgsacprofile(:,:,ii) =  avgsacprofile(:,:,ii) - min(avgsacprofile(:,:,ii));
%         avgsacprofile(:,:,ii) = avgsacprofile(:,:,ii)/max(avgsacprofile(:,:,ii));
%     end
%     
%     figure(graphnum+11)
%     subtitle('Average-Smoothed Fixation Profile by Parameter')
%     hax(11,i) = subplot(2,2,i);
%     title(tags{i})
%     hold all
%     h = area(5:fixlen-5,ones(1,fixlen-9));
%     set(h,'FaceColor',[.75 .75 .75])
%     set(h,'EdgeColor','none')
%     for ii =  1:size(avgfixprofile,3);
%         plot(avgfixprofile(:,:,ii),'linewidth',2)
%     end
%     hold off
%     xlim([1 fixlen])
%     set(gca,'XTick',[])
%     set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
%     if i == 1
%         ylabel('Normalized Value')
%         set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
%     end
%     if i == 2
%         legend([{'fixation'} variables],'Location','NorthEastOutside');
%     end
%     if i > 2
%         xlabel('Warped Time')
%     end
%     
%     figure(graphnum+12)
%     subtitle('Average-Smoothed Saccade Profile by Parameter')
%     hax(12,i) = subplot(2,2,i);
%     title(tags{i})
%     hold all
%     h1 = area(1:5,ones(1,5));
%     set(h1,'FaceColor',[.75 .75 .75])
%     set(h1,'EdgeColor','none')
%     h2 = area(saclen-4:saclen,ones(1,5));
%     set(h2,'FaceColor',[.75 .75 .75])
%     set(h2,'EdgeColor','none')
%     for ii = 1:size(avgsacprofile,3)
%         p(ii) = plot(avgsacprofile(:,:,ii),'linewidth',2);
%     end
%     hold off
%     xlim([1 saclen])
%     set(gca,'XTick',[])
%     set(gca,'YTick',[0 1],'YTickLabel',{'0','1'})
%     if i == 1
%         ylabel('Normalized Value')
%         set(get(gca,'YLabel'),'Position',[-2 -0.2 0])
%     end
%     if i == 2
%         legend([h1 p],[{'fixation'} variables],'Location','NorthEastOutside');
%     end
%     if i > 2
%         xlabel('Warped Time')
%     end
%     
%     figure(graphnum+13)
%     subtitle('Probability of Saccade Angle Changeing > 45 Degrees')
%     hax(13,i) = subplot(2,2,i);
%     title(tags{i})
%     hold on
%     h  = area(allview{i}.mediansac+1:size(allview{i}.persistence,2),...
%         ones(1,size(allview{i}.persistence,2)-allview{i}.mediansac));
%     set(h,'FaceColor',[.75 .75 .75])
%     set(h,'EdgeColor','none')
%     p = plot(mean(allview{i}.persistence));
%     hold off
%     xlim([1 size(allview{i}.persistence,2)])
%     if i == 1
%         ylabel('Probability of Saccade Angle Changeing > 45 Degrees')
%         set(get(gca,'YLabel'),'Position',[-3 -0.3 0])
%     end
%     if i == 2
%         legend([h p],{'fixation','persistence'},'Location','NorthEastOutside');
%     end
%     if i > 2
%         xlabel('Warped Time')
%     end
%     
%     figure(graphnum+14)
%     subtitle('Saccade and Fixation Distance')
%     hax(14,i) = subplot(2,2,i);
%     plot(nanmean(allview{i}.distanceprofile))
%     title(tags{i})
%     xlabel('Warped Time')
%     ylabel('Distance (pixels)')
%     
%     [allprobangle2fix] = hist(allview{i}.sacangle_2fix(~isnan(allview{i}.sacangle_2fix)),360);
%     allprobangle2fix = [allprobangle2fix(36:-1:1) allprobangle2fix allprobangle2fix(end:-1:end-36)];
%     allprobangle2fix = filtfilt(1/6*ones(1,6),1,allprobangle2fix);
%     allprobangle2fix = allprobangle2fix(37:end-37);
%     allprobangle2fix = allprobangle2fix/sum(allprobangle2fix);
%     allprobangle2fix = [allprobangle2fix allprobangle2fix(1)];
%     
%     figure(graphnum+15)
%     subtitle('Distribution of saccade angles entering fixations')
%     hax(15,i) = subplot(2,2,i);
%     polar(n,allprobangle2fix)
%     title(tags{i})
%     ph=findall(gca,'type','text');
%     set(ph,'fontweight','bold');
%     
%     figure(graphnum+16)
%     subtitle('Distribution of Saccade Amplitudes')
%     hax(16,i) = subplot(2,2,i);
%     sacdist = allview{i}.sacamplitude(~isnan(allview{i}.sacamplitude));
%     sacdist(sacdist > 800) = [];
%     hist(sacdist,100)
%     xlabel('Distance (Pixels)')
%     title(tags{i})
% end
% 
% screen_size = get(0, 'ScreenSize');
% figdir = ['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Figures\BCRW_behavior\'];
% figuretitles = {
%     'Distribution of Angles between fixations';
%     'Distribution of Angles Leaving a fixation';
%     'Distribution of Fixation Durations';
%     'Distribution of Distances Between Fixations';
%     'Fixation (Saccade) Rate';
%     'Distribution of Saccade Durations';
%     'Distribution of Saccade Distances';
%     'Correlation between Saccade Duration and Saccade Distance';
%     'Fixation and Saccade Statistics by Fixation or Saccade Number';
%     '2-D Fixation PDF'
%     'Average-Smoothed Fixation Profile';
%     'Average-Smoothed Saccade Profile';
%     'Persistence Profile';
%     'Saccade and Fixation Distance Profiles';
%     'Distribution of angles entering a fixation';
%     'Distribution of Saccade Amplitudes';
%     };
% 
% for ff = 1:size(hax,1);
%     figure(graphnum+ff)
%     set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%     for i = 1:4
%         switch ff
%             case {10,11,12,13}
%                 if i == 1;
%                     pos = [.13 0.57 0.33 0.33];
%                 elseif i == 2;
%                     pos = [.5 0.57 0.33 0.33];
%                 elseif i == 3;
%                     pos = [.13 0.13 0.33 0.33];
%                 else
%                     pos = [.5 0.13 0.33 0.33];
%                 end
%             otherwise
%                 if i == 1;
%                     pos = [.13 0.57 0.33 0.33];
%                 elseif i == 2;
%                     pos = [.57 0.57 0.33 0.33];
%                 elseif i == 3;
%                     pos = [.13 0.13 0.33 0.33];
%                 else
%                     pos = [.57 0.13 0.33 0.33];
%                 end
%         end
%         set(hax(ff,i),'position',pos);
%     end
%     h = findobj(gcf,'Type','axes','Tag','legend');
%     if ~isempty(h)
%         legpos = get(h,'Position');
%         legpos(1) = 0.875;
%         legpos(2) = 0.75;
%         set(h,'Position',legpos)
%     end
%     saveas(gcf,[figdir figuretitles{ff}])
%     print(gcf,'-r300','-djpeg',[figdir figuretitles{ff}])
% end
% 
% allviewvariables = {
%     'tags: subject names';
%     '';
%     'allview: all combined data by subject';
%     'allview.densitymap: positions of fixations';
%     'allview.allfixations: fixation profile by parameters warped twice to median length';
%     'allview.allsaccades: saccade profile by parameters warped twice to median length';
%     'allview.persistence: persistence/probability of eye movement changing >45 angle double warped';
%     'allview.anglebtwfix: angles between fixations by fixation';
%     'allview.sacangle_2fix: saccade angles entering a fixation';
%     'allview.distanceprofile: velocity of eye movements for saccade + subsequent fixation';
%     'allview.distbtwnfix: distance between fixations by fixation number';
%     'allview.fixduration: fixation duration by fixation number';
%     'allview.sacangle: angle of saccade leaving a fixation by saccade number';
%     'allview.sacdist: distance of saccade by saccade number';
%     'allview.sacamplitude: saccade amplitude from end to end';
%     'allview.timebtwfix: time between fixations by fixation number';
%     };
% 
% save(['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\'...
%     'Image Sets\CombinedViewingBehavior.mat'],'tags','allview',...
%     'allStatsbyfixation','allviewvariables');
%%
%---[4] Run simulations for manipulated images---%
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\';
% tags = {'PW','TT','TO','RR'};
% 
% Combinedbehaviorfile = ['C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\'...
%     'Image Sets\CombinedViewingBehavior.mat'];
% load(Combinedbehaviorfile,'allview')
% imageX = 800; imageY = 600;
% plotoptions.runs = 'none'; %all/none
% plotoptions.probdens = 'none';
% plotoptions.type = 'sal'; %sal/image
% IOR_tau = 1/17;
% 
% for SET =1:40;
%     if SET < 10
%         setzero = '0';
%     else
%         setzero = '';
%     end
%     
%     cd([scm_image_dir  'SCM' setzero (num2str(SET)) '\']);
%     
%     matfiles = what;
%     statfiles = [];
%     for i = 1:length(matfiles.mat);
%         str = strfind(matfiles.mat{i},'saliencemap');
%         if ~isempty(str)
%             for t = 1:length(tags)
%                 disp(['Running ' tags{t} ' on image #' num2str(i) ' from ' num2str(SET)])
%                 run_BCRWCF(allview{t},matfiles.mat{i},tags{t},imageX,imageY,plotoptions,IOR_tau)
%             end
%         end
%     end
% end
%%
%---[5] Determine the number of Predicted fixations to Manipulated ROI---%
% scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\';
% scm_eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
% scm_ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\';
% tags = {'PW','TT','TO','RR'};
% 
% type = ['rm']; %object replaced and object moved
% 
% replaced_difficulty = NaN(40,36);
% moved_difficulty = NaN(40,36);
% 
% replaced_total = NaN(40,36);
% moved_total = NaN(40,36);
% replaced_area = NaN(40,36);
% moved_area = NaN(40,36);
% %area should be the same since the same object is being manipulated

% for SET =1:40;
%     if SET < 10
%         setzero = '0';
%     else
%         setzero = '';
%     end   
%     cd([scm_image_dir  'SCM' setzero (num2str(SET)) '\']);
%     
%     load([scm_ROI_dir 'SCMRM' setzero' num2str(SET) '_ROIs.mat']);
%     
%     
%     for index = 1:36
%         if index < 10
%             indexzero = '0';
%         else
%             indexzero = '';
%         end
%         if ~exist([tags{1} '-S' setzero num2str(SET) 'I' indexzero  num2str(index) 'r-BCRW.mat'],'file')...
%                 && ~exist([tags{1} '-S' setzero num2str(SET) 'I' indexzero  num2str(index) 'm-BCRW.mat'],'file');
%             %if it doesn't exist try the next image
%             continue
%         elseif exist([tags{1} '-S' setzero num2str(SET) 'I' indexzero  num2str(index) 'r-BCRW.mat'],'file')
%             s = 1;
%         else exist([tags{1} '-S' setzero num2str(SET) 'I' indexzero  num2str(index) 'm-BCRW.mat'],'file')
%             s = 2;
%         end
%         
%         fixations = cell(1,length(tags));
%         for t = 1:length(tags)
%             load([tags{t} '-S' setzero num2str(SET) 'I' indexzero  num2str(index) type(s) '-BCRW.mat'],'fixationtimes');
%             for i = 1:size(fixationtimes,1);
%                 fixations{t}{i} = NaN(2,60);
%             end
%             for i = 1:size(fixationtimes,1);
%                 tind = find(fixationtimes(i,:,1) > 0);
%                 if length(tind) > 60
%                     tind = tind(1:60);
%                 end
%                 for ii = 1:length(tind)
%                     x = fixationtimes(i,tind(ii),1);
%                     y = fixationtimes(i,tind(ii),2);
%                     fixations{t}{i}(:,ii) = [x;y];
%                 end
%             end
%         end
%         
%         difficulty = NaN(4,length(fixations{1}));
%         total = NaN(4,length(fixations{1}));
%         for t = 1:length(fixations)
%             for i = 1:length(fixations{t})
%                 %add 1 dva or 24 pixel buffer around ROI
%                 if s==1
%                     ind = find(fixations{t}{i}(1,:)+24 > ROIs{index}{1}(1)...
%                         & fixations{t}{i}(1,:)-24 < ROIs{index}{1}(2) ....
%                         & fixations{t}{i}(2,:)+24 > ROIs{index}{1}(3) ...
%                         & fixations{t}{i}(2,:)-24 < ROIs{index}{1}(4));%BCRW y coordinates already flipped
%                 else
%                     ind = find(fixations{t}{i}(1,:)+24 >  ROIs{index}{2}(1)...
%                         & fixations{t}{i}(1,:)-24 < ROIs{index}{2}(2) ....
%                         & fixations{t}{i}(2,:)+24 > ROIs{index}{2}(3) ...
%                         & fixations{t}{i}(2,:)-24 < ROIs{index}{2}(4));%BCRW y coordinates already flipped
%                 end
%                 if ~isempty(ind)
%                     difficulty(t,i) = min(ind);
%                     total(t,i) = length(ind);
%                 else
%                     difficulty(t,i) = NaN; %if not found cap
%                     total(t,i) = 0;
%                 end
%             end
%         end
%         img = zeros(600,800);
%         for t = 1:length(fixations)
%             for i = 1:length(fixations{t})
%                 for ii = 1:size(fixations{t}{i},2)
%                     if ~isnan(fixations{t}{i}(1,ii))
%                         img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii)) =  img(fixations{t}{i}(2,ii),fixations{t}{i}(1,ii))+1;
%                     end
%                 end
%             end
%         end
%         if s == 1
%             replaced_difficulty(SET,index) = nanmean(difficulty(1:end));
%             replaced_total(SET,index) = nanmean(total(1:end));
%             replaced_area(SET,index) = (ROIs{index}{1}(2)-ROIs{index}{1}(1))*...
%                 (ROIs{index}{1}(4)-ROIs{index}{1}(3));
%         else
%             moved_difficulty(SET,index) = nanmean(difficulty(1:end));
%             moved_total(SET,index) = nanmean(total(1:end));
%             moved_area(SET,index) = (ROIs{index}{2}(2)-ROIs{index}{2}(1))*...
%                 (ROIs{index}{2}(4)-ROIs{index}{2}(3));
%         end
%     end
% end
% 
% nanmean(nanmean(replaced_difficulty))
% nanmean(nanmean(moved_difficulty))
% nanmean(nanmean(replaced_total))
% nanmean(nanmean(moved_total))
% nanmean(nanmean(replaced_area))
% nanmean(nanmean(moved_area))
% 
% [~,tpd] = ttest2(replaced_difficulty(1:end),moved_difficulty(1:end))
% [~,kpd] = kstest2(replaced_difficulty(1:end),moved_difficulty(1:end))
% [~,tpt] = ttest2(replaced_total(1:end),moved_total(1:end))
% [~,kpt] = kstest2(replaced_total(1:end),moved_total(1:end))
% [~,tpa] = ttest2(replaced_area(1:end),moved_area(1:end))
% [~,kpa] = kstest2(replaced_area(1:end),moved_area(1:end))

%%
%---[6] Determine the number of Observed fixations to Manipulated ROI---%
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\';
scm_eyedata_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Eye Data\';
scm_ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\';
tags = {'PW','TT','TO','RR'};

type = ['rm']; %object replaced and object moved


moved_fixations = NaN(40,36,length(tags)); %proportion total number of fixations
moved_fix10 = NaN(40,36,length(tags)); %proportion of of 1st 10 fixations
moved_time = NaN(40,36,length(tags)); %proportion of total time
moved_time3 = NaN(40,36,length(tags));%proprotion of 1st 3 seconds
moved_time15 = NaN(40,36,length(tags));%proprotion of 1st 3 seconds

replaced_fixations = NaN(40,36,length(tags));  %proportion total number of fixations
replaced_fix10 = NaN(40,36,length(tags));%proportion of of 1st 10 fixations
replaced_time = NaN(40,36,length(tags));%proportion of total time
replaced_time3 = NaN(40,36,length(tags));%proprotion of 1st 3 seconds

cd(scm_eyedata_dir);
matfiles = what;
eyedatafiles = [];
for i = 1:length(matfiles.mat);
    str = strfind(matfiles.mat{i},'ROIdata.mat');
    if ~isempty(str)
        monk = [];
        for t = 1:4
            if strcmpi(tags{t},matfiles.mat{i}(1:2))
                monk = t;
                break
            end
        end
        
        load(matfiles.mat{i})
        set = str2double(itemfile(6:7));
        
        for im = 1:length(image_numbers{3}) %replaced
            replaced_fixations(set,image_numbers{3}(im),monk) = allnumFixationsInROI{3}(im);%proportion total number of fixations
            replaced_fix10(set,image_numbers{3}(im),monk) = num10FixationsInROI{3}(im);%proportion of of 1st 10 fixations
            replaced_time(set,image_numbers{3}(im),monk) = allTime_ROI_timeWindow{3}(im);%proportion of total time
            replaced_time3(set,image_numbers{3}(im),monk) = Time3_ROI_timeWindow{3}(im);%proprotion of 1st 3 seconds
        end
        
        for im = 1:length(image_numbers{4}) %moved
            moved_fixations(set,image_numbers{4}(im),monk) = allnumFixationsInROI{4}(im);%proportion total number of fixations
            moved_fix10(set,image_numbers{4}(im),monk) = num10FixationsInROI{4}(im);%proportion of of 1st 10 fixations
            moved_time(set,image_numbers{4}(im),monk) = allTime_ROI_timeWindow{4}(im);%proportion of total time
            moved_time3(set,image_numbers{4}(im),monk) = Time3_ROI_timeWindow{4}(im);%proprotion of 1st 3 seconds
        end
        
    end
end
%%

moved_fixations = nanmean(moved_fixations,3);
moved_fix10 = nanmean(moved_fix10,3);
moved_time = nanmean(moved_time,3);
moved_time3 = nanmean(moved_time3,3);
moved_time15 = nanmean(moved_time15,3);

replaced_fixations = nanmean(replaced_fixations,3);
replaced_fix10 = nanmean(replaced_fix10,3);
replaced_time = nanmean(replaced_time,3);
replaced_time3 = nanmean(replaced_time3,3);


temp1 = moved_fixations(1:end);
temp2 = moved_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_fixations(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end


temp1 = moved_fix10(1:end);
temp2 = moved_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2); 

figure
plot(moved_fix10(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time(1:end);
temp2 = moved_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2); 
figure
plot(moved_time(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time3(1:end);
temp2 = moved_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2); 
figure
plot(moved_time3(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time3(1:end);
temp2 = moved_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2); 
figure
plot(moved_time15(1:end),moved_difficulty(1:end),'k.')
xlabel('Proportion of 0.5-1.5 secs')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end


temp1 = moved_fixations(1:end);
temp2 = moved_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_fixations(1:end),moved_total(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_fix10(1:end);
temp2 = moved_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_fix10(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time(1:end);
temp2 = moved_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_time(1:end),moved_total(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time3(1:end);
temp2 = moved_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_time3(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end

temp1 = moved_time15(1:end);
temp2 = moved_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(moved_time15(1:end),moved_total(1:end),'k.')
xlabel('Proportion of 0.5-1.5 secs')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Moved r = ' num2str(r(2))]);
else
    title('Moved')
end


temp1 = replaced_fixations(1:end);
temp2 = replaced_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_fixations(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_fix10(1:end);
temp2 = replaced_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_fix10(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_time(1:end);
temp2 = replaced_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_time(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_time3(1:end);
temp2 = replaced_difficulty(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_time3(1:end),replaced_difficulty(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('First Fixation in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_fixations(1:end);
temp2 = replaced_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_fixations(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of All Fixations')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_fix10(1:end);
temp2 = replaced_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_fix10(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of 1st 10 Fixations')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end

temp1 = replaced_time(1:end);
temp2 = replaced_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_time(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of All Time')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end


temp1 = replaced_time3(1:end);
temp2 = replaced_total(1:end);
temp1(isnan(temp2)) = [];
temp2(isnan(temp2)) = [];
temp2(isnan(temp1)) = [];
temp1(isnan(temp1)) = [];
[r,p] = corrcoef(temp1,temp2);

figure
plot(replaced_time3(1:end),replaced_total(1:end),'k.')
xlabel('Proportion of 1st 3 secs')
ylabel('Total Fixations in ROI from BCRW')
if p(2) < 0.05/10
    title(['Replaced r = ' num2str(r(2))]);
else
    title('Replaced')
end
