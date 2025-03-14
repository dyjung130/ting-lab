%This is a script for analyzing Bolt et al., 2024 data set
%(ME-...something) 030825 DJ
% For the EEG electrode location, using the montage from
% brainproducts-RNP-BA-128: Brain Products with 10-10 electrode names (128 channels) 
% may work.
%%
clear all
%% check EEG montage
dirlist = "/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/matlab-lib/fieldtrip-20250114/template/electrode/";
filename = fullfile(dirlist,'standard_1020.elc');
standard_1020 = ft_read_sens(filename);
%% setup directories
baseDir = '/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/bolt2025/dataset_chang_cue';%/eeg/raw';
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


%% iterate for each EEG file
[E_ft,E_ft_laplace] = deal(cell(length(eegFiles),1));

for i = 1:length(eegFiles)
    E = load(fullfile(baseDir,eegDir,eegFiles{i}));% I think this is sEEGLAB format.
    %the loaded data have the channel location info in E.EEG.chanlocs.
    addpath('lib/');%add library
    % Convert to FieldTrip format for analysis
    Etemp = struct();
    Etemp.fsample = E.EEG.srate;%sampling rate (@ 250Hz)
    Etemp.hdr = struct();
    Etemp.trial = {E.EEG.data};%data
    Etemp.time = {E.EEG.times./1000};%time (the time in the original data appears to be in ms)
    Etemp.label = {E.EEG.chanlocs.labels}';% channel labels
    Etemp.hdr.Fs = E.EEG.srate;%set to (@ 250Hz)
    Etemp.hdr.nChans = length(Etemp.label);%number of channels
    Etemp.hdr.label = Etemp.label;%
    Etemp.cfg = struct();
    E_ft{i} = Etemp;
    % Using the standard 10-20 template, provided from MNE (FieldTrip), which has the somewhat acutal phyiscal coordinates*
    %the standard_1020 file is in 'mm'
    eegLocInd = cellfun(@(x) (cellfun(@(y) strcmpi(x,y),Etemp.label,'uni',0)), standard_1020.label,'uni',0);
    eegLocInd = find(sum(cell2mat(horzcat(eegLocInd{:}))));
    standard_1020_mod.chanpos = standard_1020.chanpos(eegLocInd,:);
    standard_1020_mod.chantype = standard_1020.chantype{eegLocInd};
    standard_1020_mod.chanunit = standard_1020.chanunit(eegLocInd);
    standard_1020_mod.elecpos = standard_1020.elecpos(eegLocInd,:);
    standard_1020_mod.label = standard_1020.label(eegLocInd);%this will narrow down to the 32 channel 10-20 system.
    % Laplace transform (sensor space)
    %for fieldtrip setup
    actualElecPos = [];
    actualElecPos.label = standard_1020_mod.label;
    actualElecPos.pnt = standard_1020_mod.chanpos;
    %
    cfg = [];
    cfg.layout.pos = actualElecPos.pnt(:,[1,2]);% X and Y positions
    cfg.layout.label = actualElecPos.label;
    lay = ft_prepare_layout(cfg);
    cfg_neighb = [];
    cfg_neighb.method = 'triangulation';
    cfg_neighb.layout = lay;
    neighbours  = ft_prepare_neighbours(cfg_neighb, Etemp);

    % laplace transform it
    cfg = [];
    cfg.method = 'spline';%'finite';
    cfg.degree = 7;
    cfg.elec         = actualElecPos;%structure with electrode positions or filename, see FT_READ_SENS
    cfg.trials       = 'all';% or a selection given as a 1xN vector (default = 'all')
    cfg.neighbours   = neighbours; %neighbourhood structure, see FT_PREPARE_NEIGHBOURS
    [E_ft_laplace{i}] = ft_scalpcurrentdensity(cfg, Etemp);
end
generatePowerSpectra(E_ft_laplace,[1,100],[1,100]);
return;
%% Topographical illustration
%display the top-view coordinate of the EEG
figure;
theta = linspace(0, 2*pi, 100); x = cos(theta); y = sin(theta); plot(x,y,'k','linewidth',1); hold on;
subplot(1,3,1);
scatter(actualElecPos.pnt(:,1),actualElecPos.pnt(:,2),'or','filled');
text(actualElecPos.pnt(:,1),actualElecPos.pnt(:,2)+0.1,actualElecPos.label);
%xlim([-1.25 1.25]); ylim([-1.25,1.25]);
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
E_ft{i}.avg{1} = E_ft{i}.trial{1};
ft_topoplotER(cfg,E_ft{i});
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
E_ft_laplace{i}.avg{1} = E_ft_laplace{i}.trial{1};
ft_topoplotER(cfg,E_ft_laplace{i});
pbaspect([1 1 1]);
title('EEG, Sensor-space');

maxClim = max(abs(caxis));
caxis([-maxClim,maxClim]);

%match_clim(ax);

%% MRI source localization
% %first find the subject from the anatomical result (find index)
matched = find(cellfun(@(x) length(regexp(x,eegSubjects{i})), anatSubjects));%index of the anat mri file

%check if there is only one file
if length(matched) ~= 1
    warning('There are more than 1 file with the same name.');
    matched = matched(1); %since it will be identical file (hopefully), just grab one of the file
end

%load the headmodel from coregisterMRIHeadmodelWithEEGSensors.m function
dirHeadmodels = dir('mri');
allHeadmodels =  {dirHeadmodels.name};
headmodel = cellfun(@(x) regexp(x,strcat('.*',anatSubjects{matched},'.*'),'match'),allHeadmodels,'uni',0);
headmodel = headmodel{~cellfun(@isempty,headmodel)};
%headmodel = load(fullfile('mri',headmodel{1}));
headmodel = load(fullfile('mri','mriModel_51728393_sub_0020-mprage'));
headmodel = headmodel.vol;

cfg = [];
cfg.elec = standard_1020_mod;
cfg.headmodel = headmodel;
%cfg.reducerank = 3;
cfg.resolution = 10;
cfg.sourcemodel.unit = 'mm';
sourcemodel = ft_prepare_leadfield(cfg);
%%
set(0,'DefaultFigureWindowStyle','docked')
%shorten data first (for testing)
cfg = [];
cfg.latency = [600 700];%in seconds
cfg.channel = ft_channelselection('eeg',E_ft{i});% just grab the EEG channels (no ECG..!)
E_temp = ft_selectdata(cfg,Etemp);

%re-reference to the average of the channels
cfg = [];
cfg.channel = 'all';
cfg.reref = 'yes';
cfg.refmethod = 'avg';%common 
cfg.refchannel = 'all';
E_temp = ft_preprocessing(cfg,E_temp);

%for DICS: Dynamical Imaging of Coherent Sources (DICS) (Gross et al. 2001)
%for (oscillatory) power
% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.output    = 'powandcsd';
% cfg.tapsmofrq = 4;
% cfg.foilim    = [10 10];
% freqPost = ft_freqanalysis(cfg, E_temp);
% freqPost.elec = standard_1020_mod;

% %for LCMV: Linearly Constrained Minimum Variance (LCMV) and the estimates are calculated in the time domain (van Veen et al., 1997).
% %for amplitude
cfg = [];
cfg.covariance = 'yes';
E_timelock = ft_timelockanalysis(cfg,E_temp);
% %for quick plot
% %figure, plot(E_timelock.time, E_timelock.avg);

%%
% cfg              = [];
% %cfg.method       = 'lcmv';
% cfg.method = 'dics';
% cfg.frequency    = 10;%for DICS
% cfg.sourcemodel  = sourcemodel;
% cfg.headmodel    = headmodel;
% cfg.dics.projectnoise = 'yes';
% cfg.dics.lambda = 1;
% %cfg.lcmv.keepfilter  = 'yes';
% %cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
% sourcePost_nocon = ft_sourceanalysis(cfg, freqPost);
%https://eeglab.org/tutorials/09_source/EEG_sources.html
cfg              = [];
cfg.method       = 'lcmv';
cfg.sourcemodel  = sourcemodel;
cfg.headmodel    = headmodel;
source = ft_sourceanalysis(cfg, E_timelock);

cfg = [];
cfg.projectmom = 'yes';
cfg.flipori = 'yes';
sourceProj = ft_sourcedescriptives(cfg, source);

cfg = [];
cfg.parameter = 'mom';%dipole moment
cfg.operation = 'abs';
sourceProj = ft_math(cfg, sourceProj);

cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'mom';
figure; ft_sourceplot(cfg, sourceProj);
%%
fileToRead = fullfile(dirName,fileName);
mri = ft_read_mri(fileToRead,'dataformat','nifti');
mri = ft_convert_coordsys(mri,'mni','feedback','no');
%%
cfg              = [];
cfg.downsample   = 1;
cfg.parameter    = 'pow';
source.oridimord = 'pos';
source.momdimord = 'pos';
sourceInt  = ft_sourceinterpolate(cfg, source , mri);

cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
ft_sourceplot(cfg, sourceInt);
%%
% cfg            = [];
% cfg.downsample = 1;
% cfg.parameter  = 'pow';
% sourcePostInt_nocon  = ft_sourceinterpolate(cfg, source, mri);

%%
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
%cfg.funcolorlim = [0.0 1.2];
%cfg.opacitymap = [0.0 1.2];
cfg.opacitymap = 'rampup'
cfg.maskparameter = cfg.funparameter;
ft_sourceplot(cfg, sourceInt);
%% FOOOF analysis on a single representative channel

% 
% set(0,'DefaultFigureWindowStyle','docked')
% targetChannelIndex = 20;
% %shorten data first (for testing)
% cfg = [];
% cfg.latency = [1 100];%in seconds
% cfg.channel = targetChannelIndex;%channe of interest
% E_zscored_shorten = ft_selectdata(cfg,E_ft_laplace{i});
% 
% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.output    = 'fooof_aperiodic';
% cfg.pad='nextpow2';
% cfg.tapsmofrq = 4;
% cfg.foilim    = [1 100]; %Hz
% E_1overf = ft_freqanalysis(cfg, E_zscored_shorten);%aperiodic component (fractal)
% 
% cfg.output = 'pow';
% E_orig = ft_freqanalysis(cfg,E_zscored_shorten);%original spectrum
% 
% %subtract the 1/f component from the power spectrum
% cfg =  [];
% cfg.parameter = 'powspctrm';
% cfg.operation = 'x2./x1';
% E_osc = ft_math(cfg,E_1overf,E_orig);%oscillations
% 
% figure(2);
% plot(log(E_orig.freq),log(E_orig.powspctrm),'b','LineWidth',2);
% hold on;
% plot(log(E_1overf.freq),log(E_1overf.powspctrm),'k','LineWidth',2);
% plot(log(E_osc.freq),log(E_osc.powspctrm),'r','LineWidth',2)
% xlabel('Frequency (log(Hz))');
% ylabel('Power (log(arb. unit))');
% %draw some vertical lines for better understanding of the plot
% fRange = [4,8,12,30,60,100];
% yline(1,':'); %boundary for the x2/x1 (for actual oscillatory component
% xline(log(fRange),'--');
% text(log(fRange),ones(1,length(fRange))*3,cellstr(strtrim(num2str(fRange'))));
% xlim([0,ceil(log(max(fRange)))]);
% ylim([min(ylim),max(ylim)+1.1]);
% title(E_zscored_shorten.label{1});
% pbaspect([1 1 1]);
% legend({'Original','Aperiodic','Oscillatory'},'Location','SouthWest');
% hold off;
% 
