% Code Written by Seth Koenig 9/1/2015
% code automatically runs all analysis for SCMRM the Scene Manipulation
% task for the relational memory project
clear,clc


cortex_files = {'TO160629.2','TO160630.2','TO160705.2'};
item_files = {'SCMRM17.itm','SCMRM18.itm','SCMRM19.itm'};
heads_free = [1 1 0];


fixwin = 3.5; %size of the fixation window/2. Width of fixwin is 7 dva
imageX = 800;
imageY = 600;

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\\Tobii SRI eye tracking\Eye Data\';%location of processed eye data
ROI_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\'; %locations of regions of interest
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Tobii SRI eye tracking\Figures\'; %where to put plot
img_dir = 'C:\Users\seth.koenig\Documents\MATLAB\SCM_RelationalMemory\Image Sets\'; %where all the image sets are


%%---Preprocess all the ListRM data---%%%
for file = 1:length(cortex_files)
    ImportSCMRMData(cortex_files{file},item_files{file},figure_dir,data_dir);
end

%%---Determine Amount of Time/# of Fixations in ROI---%%%
for file = 1:length(cortex_files)
    SCMRM_ROIanalysis(cortex_files{file},item_files{file},figure_dir,data_dir,ROI_dir,img_dir);
end
%%
%%%---Combined Data across sets---%%%
CombineSCMRM_ROIdata(cortex_files,data_dir,figure_dir)
%%
%%%---Determine Fixation Durations, Saccade Amplitudes, Pupil Diameter Across Multiple Sets---%
combinedEyeMovementStats_SCMRM(data_dir,cortex_files)

% Pupil correlated with behavior
combinedPupilROI_SCMRM(data_dir,cortex_files,ROI_dir)