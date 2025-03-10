baseDir = '/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/bolt2025/dataset_chang';%/eeg/raw';
eegDir = '/eeg/raw/';% EEG
anatDir = '/anat/raw/';%anatomical MRI
funcDir = '/func/raw/';%functional MRI
%% Just sort data based on file names
anatFiles = dir(fullfile(baseDir,anatDir));
anatFiles = {anatFiles.name};
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*_sub_*','match')), anatFiles, 'uni', 0));
anatFiles = anatFiles(logical(fileIndex));
anatSubjects = cellfun(@(x) regexp(x,'sub_[0-9]*','match'),anatFiles);% available subjects identified from anat MRI
%% Now sort EEG data based on the file names obtained from anatomical...yeee
eegFiles = dir(fullfile(baseDir,eegDir));
eegFiles = {eegFiles.name};
%divide into "subject + mr session" 
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*sub_[0-9]*-*')), eegFiles, 'uni', 0));
eegFiles = eegFiles(logical(fileIndex));
eegSubjects = cellfun(@(x) regexp(x,'sub_[0-9]*','match'),eegFiles);%available subjects for EEG (can be more files than the identified number of subjects)
%%
i = 1;%eeg file indexs
s
E = load(fullfile(baseDir,eegFiles{i}));% I think this is sEEGLAB format.
%the loaded data have the channel location info in E.EEG.chanlocs.
addpath('lib/');%add library
%% Convert to FieldTrip format for analysis
E_ft = struct();
E_ft.fsample = E.EEG.srate;%sampling rate (@ 250Hz)
E_ft.hdr = struct();
E_ft.trial = {E.EEG.data};%data 
E_ft.time = {E.EEG.times./1000};%time (the time in the original data appears to be in ms)
E_ft.label = {E.EEG.chanlocs.labels}';% channel labels
E_ft.hdr.Fs = E.EEG.srate;%set to (@ 250Hz)
E_ft.hdr.nChans = length(E_ft.label);%number of channels
E_ft.hdr.label = E_ft.label;%
E_ft.cfg = struct();

%% Laplace transform (sensor space)
%POLHEMUS_ELEC_FILE = 'C:\Users\dyjun\Downloads\w-tdoaprob03\EEG_coordinate\whopper_20211103.txt';
actualElecPos.label = E_ft.label;%channel labels
%since last channel is ECG (electrocardiogram), remove it
actualElecPos.label(end) = [];
actualElecPos.pnt = vertcat([E.EEG.chanlocs.X], [E.EEG.chanlocs.Y], [E.EEG.chanlocs.Z])';%transposed
%%Since the montage/layout is from Posterior to Anterior (left to right), rotate -90
%degrees so that the montage/layout matches the direction of the face and
%ears in the output figure 

% Create rotation matrix
theta = -90; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
actualElecPos.pnt(:,[1,2]) = actualElecPos.pnt(:,[1,2]) * R;



%%
cfg = [];
cfg.layout.pos = actualElecPos.pnt(:,[1,2]);% X and Y positions
cfg.layout.label = actualElecPos.label;
lay = ft_prepare_layout(cfg);
cfg_neighb = [];
cfg_neighb.method = 'triangulation';
cfg_neighb.layout = lay;
neighbours  = ft_prepare_neighbours(cfg_neighb, E_ft);

% laplace transform it
cfg = [];
cfg.method = 'spline';%'finite';
cfg.degree = 7;
cfg.elec         = actualElecPos;%structure with electrode positions or filename, see FT_READ_SENS
cfg.trials       = 'all';% or a selection given as a 1xN vector (default = 'all')
cfg.neighbours   = neighbours; %neighbourhood structure, see FT_PREPARE_NEIGHBOURS
[E_ft_laplace] = ft_scalpcurrentdensity(cfg, E_ft);
%% Topographical illustration
%display the top-view coordinate of the EEG
figure(1)
theta = linspace(0, 2*pi, 100); x = cos(theta); y = sin(theta); plot(x,y,'k','linewidth',1); hold on;
subplot(1,3,1);
scatter(actualElecPos.pnt(:,1),actualElecPos.pnt(:,2),'or','filled');
text(actualElecPos.pnt(:,1),actualElecPos.pnt(:,2)+0.1,E_ft.label(1:31));
xlim([-1.25 1.25]); ylim([-1.25,1.25]);
pbaspect([1 1 1]);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('EEG Montage');

ax = [];

%plot raw (processed) eeg
ax(1) = subplot(1,3,2);
cfg = [];
cfg.layout = lay;
cfg.trials = 1;
cfg.xlim = [0 0];
cfg.interpolation = 'V4';
cfg.interactive = 'no';
cfg.colorbar = 'EastOutside';
cfg.comment = 'no';
cfg.gridscale = 200;
cfg.figure = 1;%add to figure (1)
E_ft.avg{1} = E_ft.trial{1};
ft_topoplotER(cfg,E_ft);
pbaspect([1 1 1]);
title('EEG');
maxClim = max(abs(caxis));
caxis([-maxClim,maxClim]);

%plot laplace tranform version
ax(2) = subplot(1,3,3);
cfg = [];
cfg.layout = lay;
cfg.trials = 1;
cfg.xlim = [0 0];
cfg.interpolation = 'V4';
cfg.interactive = 'no';
cfg.colorbar = 'EastOutside';
cfg.comment = 'no';
cfg.gridscale = 200;
cfg.figure = 1;%add to figure (1)
E_ft_laplace.avg{1} = E_ft_laplace.trial{1};
ft_topoplotER(cfg,E_ft_laplace);
pbaspect([1 1 1]);
title('EEG, Sensor-space');

maxClim = max(abs(caxis));
caxis([-maxClim,maxClim]);

%match_clim(ax);

%% MRI source localization
%first find the subject from the anatomical result (find index)
matched = find(cellfun(@(x) length(regexp(x,eegSubjects{i})), anatSubjects));%index of the anat mri file
%check if there is only one file 
if length(matched) ~= 1
    warning('There are more than 1 file with the same name.');
    matched = matched(1); %since it will be identical file (hopefully), just grab one of the file
end

fileToRead = fullfile(baseDir,anatDir,anatFiles{matched});
mri = ft_read_mri(fileToRead,'dataformat','nifti');
cfg           = [];
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);
%save segmentedmri segmentedmri
cfg = [];
cfg.method='singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

%Do you want to change the anatomical labels for the axes [Y, n]? Y
%What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
%What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
%What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s

save vol vol

vol = ft_convert_units(vol, 'cm');
figure;
ft_plot_headmodel(vol);

%% source localization on power
%https://www.fieldtriptoolbox.org/tutorial/source/beamformer/
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.pad='nextpow2';
cfg.tapsmofrq = 4;
cfg.foilim    = [18 18];
E_freq = ft_freqanalysis(cfg, E_ft_laplace);