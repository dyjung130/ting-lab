% This is a script for getting the head/source model for analyzing Bolt et al., 2024 data set
% For the EEG electrode location, using the montage from
% brainproducts-RNP-BA-128: Brain Products with 10-10 electrode names (128 channels)
% may work.
%%
clear all
addpath('lib/');%add library


%% setup directories
baseDir = '/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/bolt2025/dataset_chang_bh';%/eeg/raw';
eegDir = '/eeg/raw/';% EEG
anatDir = '/anat/raw/';%anatomical MRI
funcDir = '/func/raw/';%functional MRI
%% Just sort data based on file names
anatFiles = dir(fullfile(baseDir,anatDir));
anatFiles = {anatFiles.name};
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*_sub_*','match')), anatFiles, 'uni', 0));
anatFiles = anatFiles(logical(fileIndex));eegFiles = dir(fullfile(baseDir,eegDir));
eegFiles = {eegFiles.name};
%divide into "subject + mr session"
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*sub_[0-9]*-*')), eegFiles, 'uni', 0));
eegFiles = eegFiles(logical(fileIndex));
%% Convert to FieldTrip format for analysis
%just grab the first file (1) since all subjects should have used the same
%EEG cap.
E = load(fullfile(baseDir,eegDir,eegFiles{1}));% I think this is sEEGLAB format.
%the loaded data have the channel location info in E.EEG.chanlocs
channelLabel = {E.EEG.chanlocs.labels}';% channel labels
%% check EEG montage
%use the fieldtrip's EEG template
dirlist = "/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/matlab-lib/fieldtrip-20250114/template/electrode/";
filename = fullfile(dirlist,'standard_1020.elc');
standard_1020 = ft_read_sens(filename);
%from the standard 10-20 EEG template, select only the electrodes that are
%available in the EEG file (~32 channels).
eegLocInd = cellfun(@(x) (cellfun(@(y) strcmpi(x,y),channelLabel,'uni',0)), standard_1020.label,'uni',0);
eegLocInd = find(sum(cell2mat(horzcat(eegLocInd{:}))));
standard_1020_mod.chanpos = standard_1020.chanpos(eegLocInd,:);
standard_1020_mod.chantype = standard_1020.chantype{eegLocInd};
standard_1020_mod.chanunit = standard_1020.chanunit(eegLocInd);
standard_1020_mod.elecpos = standard_1020.elecpos(eegLocInd,:);
standard_1020_mod.label = standard_1020.label(eegLocInd);%this will narrow down to the 32 channel 10-20 system.

%% just try to plot MRI head model with the EEG snenosr
for i = 2%1:length(anatFiles)
    dirName = fullfile(baseDir,anatDir);
    fileName = anatFiles{i};
    getHeadmodel(dirName,fileName,standard_1020_mod,1);
end
