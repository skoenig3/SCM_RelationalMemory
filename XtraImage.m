%Import XtrImage Data---For reaclimation and condition to pictures not to
%be used for actual data analysis
%uses image sets from ListSQ but all novel images

clear,clc

cortex_files = {'TO170131.3'};
item_files = {'xtrimg01.itm'};

fixwin = 3.5; %size of the fixation window/2. Width of fixwin is 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Xtrimg\';%location of processed eye data
figure_dir = data_dir;
for file = 1:length(cortex_files)
    ImportSCMRMData(cortex_files{file},item_files{file},figure_dir,data_dir);
end
