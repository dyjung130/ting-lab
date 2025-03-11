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
baseDir = '/Users/dennis.jungchildmind.org/OneDrive - Child Mind Institute/bolt2025/dataset_chang_bh';%/eeg/raw';
eegDir = '/eeg/raw/';% EEG
anatDir = '/anat/raw/';%anatomical MRI
funcDir = '/func/raw/';%functional MRI
%% Just sort data based on file names
anatFiles = dir(fullfile(baseDir,anatDir));
anatFiles = {anatFiles.name};
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*_sub_*','match')), anatFiles, 'uni', 0));
anatFiles = anatFiles(logical(fileIndex));
anatSubjects = cellfun(@(x) regexp(x,'sub_[0-9]*','match'),anatFiles);% available subjects identified from anat MRI
%% just try to plot MRI head model with the EEG snenosr
for i = 1:length(anatFiles)
    dirName = fullfile(baseDir,anatDir);
    fileName = anatFiles{i};
    coregisterMRIHeadmodelWithEEGSensors(dirName,fileName,standard_1020_mod,0);
end
%% Now sort EEG data based on the file names obtained from anatomical...yeee
eegFiles = dir(fullfile(baseDir,eegDir));
eegFiles = {eegFiles.name};
%divide into "subject + mr session"
fileIndex = cell2mat(cellfun(@(x) length(regexp(x, '.*sub_[0-9]*-*')), eegFiles, 'uni', 0));
eegFiles = eegFiles(logical(fileIndex));
eegSubjects = cellfun(@(x) regexp(x,'sub_[0-9]*','match'),eegFiles);%available subjects for EEG (can be more files than the identified number of subjects)


%% iterate for each EEG file
E_ft_laplace = cell(length(eegFiles),1);

for i = 1:length(eegFiles)
    E = load(fullfile(baseDir,eegDir,eegFiles{i}));% I think this is sEEGLAB format.
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
    %% Using the standard 10-20 template, provided from MNE (FieldTrip), which has the somewhat acutal phyiscal coordinates*
    %the standard_1020 file is in 'mm'
    eegLocInd = cellfun(@(x) (cellfun(@(y) strcmpi(x,y),E_ft.label,'uni',0)), standard_1020.label,'uni',0);
    eegLocInd = find(sum(cell2mat(horzcat(eegLocInd{:}))));
    standard_1020_mod.chanpos = standard_1020.chanpos(eegLocInd,:);
    standard_1020_mod.chantype = standard_1020.chantype{eegLocInd};
    standard_1020_mod.chanunit = standard_1020.chanunit(eegLocInd);
    standard_1020_mod.elecpos = standard_1020.elecpos(eegLocInd,:);
    standard_1020_mod.label = standard_1020.label(eegLocInd);%this will narrow down to the 32 channel 10-20 system.
    %% Laplace transform (sensor space)
    %POLHEMUS_ELEC_FILE = 'C:\Users\dyjun\Downloads\w-tdoaprob03\EEG_coordinate\whopper_20211103.txt';
    actualElecPos.label = E_ft.label;%channel labels
    %since last channel is ECG (electrocardiogram), remove it
    actualElecPos.label(end) = [];
    %actualElecPos.pnt = vertcat([E.EEG.chanlocs.X], [E.EEG.chanlocs.Y],[E.EEG.chanlocs.Z])';%transposed (model, unit circle)
    actualElecPos.pnt = standard_1020_mod.chanpos;
    %%Since the montage/layout is from Posterior to Anterior (left to right), rotate -90
    %degrees so that the montage/layout matches the direction of the face and
    %ears in the output figure

    % % Create rotation matrix (only for the model, unit circle)
    % theta = -90; % to rotate 90 counterclockwise
    % R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    % actualElecPos.pnt(:,[1,2]) = actualElecPos.pnt(:,[1,2]) * R;



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
    [E_ft_laplace{i}] = ft_scalpcurrentdensity(cfg, E_ft);
end


%% Topographical illustration

%display the top-view coordinate of the EEG
figure(1)
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
E_ft_laplace{1}.avg{1} = E_ft_laplace{1}.trial{1};
ft_topoplotER(cfg,E_ft_laplace{1});
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

%for subject i = 9
fileToRead = fullfile(baseDir,anatDir,anatFiles{matched});
mri = ft_read_mri(fileToRead,'dataformat','nifti');
%rotate
mri = ft_convert_coordsys(mri,'mni','feedback','no');
pauseDuration = 1;
%Do you want to change the anatomical labels for the axes [Y, n]? Y
%What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
%What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
%What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
%Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), s(scanner origin), n(ot a landmark)? npause(1);

%segment it!
cfg           = [];
cfg.output    = 'brain';%'scalp';
segmentedmri  = ft_volumesegment(cfg, mri);
%save segmentedmri segmentedmri
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);


%save vol vol
vol = ft_convert_units(vol, 'mm');% I think the unit is already in mm but just in case it is not mm, do conversion
figure;
ft_plot_sens(standard_1020_mod);
hold on;
ft_plot_headmodel(vol);%results in mm

%% FOOOF analysis on a single representative channel
set(0,'DefaultFigureWindowStyle','docked')
targetChannelIndex = 20;
%shorten data first (for testing)
cfg = [];
cfg.latency = [1 100];%in seconds
cfg.channel = targetChannelIndex;%channe of interest
E_zscored_shorten = ft_selectdata(cfg,E_ft_laplace{i});

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fooof_aperiodic';
cfg.pad='nextpow2';
cfg.tapsmofrq = 4;
cfg.foilim    = [1 100]; %Hz
E_1overf = ft_freqanalysis(cfg, E_zscored_shorten);%aperiodic component (fractal)

cfg.output = 'pow';
E_orig = ft_freqanalysis(cfg,E_zscored_shorten);%original spectrum

%subtract the 1/f component from the power spectrum
cfg =  [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2./x1';
E_osc = ft_math(cfg,E_1overf,E_orig);%oscillations

figure(2);
plot(log(E_orig.freq),log(E_orig.powspctrm),'b','LineWidth',2);
hold on;
plot(log(E_1overf.freq),log(E_1overf.powspctrm),'k','LineWidth',2);
plot(log(E_osc.freq),log(E_osc.powspctrm),'r','LineWidth',2)
xlabel('Frequency (log(Hz))');
ylabel('Power (log(arb. unit))');
%draw some vertical lines for better understanding of the plot
fRange = [4,8,12,30,60,100];
yline(1,':'); %boundary for the x2/x1 (for actual oscillatory component
xline(log(fRange),'--');
text(log(fRange),ones(1,length(fRange))*3,cellstr(strtrim(num2str(fRange'))));
xlim([0,ceil(log(max(fRange)))]);
ylim([min(ylim),max(ylim)+1.1]);
title(E_zscored_shorten.label{1});
pbaspect([1 1 1]);
legend({'Original','Aperiodic','Oscillatory'},'Location','SouthWest');
hold off;

%% FOOOF for all channels
%normalize data first
E_zcored = [];
for i = 1:length(E_ft_laplace)
    E_zscored{i} = E_ft_laplace{i};
    E_zscored{i}.trial = cellfun(@(x) zscore(x,[],'all'),E_zscored{i}.trial,'uni',0);
end

set(0,'DefaultFigureWindowStyle','docked')
figure(3);
for targetChannelIndex = 1:length(E_ft_laplace{1}.label)
    ax_all = [];
    E_osc = [];
    for i = 1:length(E_ft_laplace)

        %shorten data first (for testing)
        cfg = [];
        cfg.latency = [1 100];%in seconds
        cfg.channel = targetChannelIndex;%channe of interest
        E_zscored_shorten = ft_selectdata(cfg,E_zscored{i});

        cfg = [];
        cfg.method    = 'mtmfft';
        cfg.output    = 'fooof_aperiodic';
        cfg.pad='nextpow2';
        cfg.tapsmofrq = 4;
        cfg.foilim    = [1 100]; %Hz
        E_1overf = ft_freqanalysis(cfg, E_zscored_shorten);%aperiodic component (fractal)

        cfg.output = 'pow';
        E_orig = ft_freqanalysis(cfg,E_zscored_shorten);%original spectrum

        %subtract the 1/f component from the power spectrum
        cfg =  [];
        cfg.parameter = 'powspctrm';
        cfg.operation = 'x2./x1';
        E_osc{i} = ft_math(cfg,E_1overf,E_orig);%oscillations
    end

    powerSpectra = cellfun(@(x) x.powspctrm, E_osc, 'uni',0);
    powerSpectra = log(vertcat(powerSpectra{:}));%concatenate (subjects x freqpow)
    %
    meanSpectrum = mean(powerSpectra,1);
    stdSpectrum = std(powerSpectra,[],1);
    semSpectrum = stdSpectrum/sqrt(size(powerSpectra,1)-1);

    ax_all(targetChannelIndex) = subplot(4, ceil(length(E_ft_laplace{i}.label)/4),targetChannelIndex);
    hold on;
    %plot(log(E_osc.freq),log(E_osc.powspctrm),'LineWidth',2)
    shadedErrorBar(log(E_osc{1}.freq),(meanSpectrum),(semSpectrum))
    %xlabel('Frequency (log(Hz))');
    %ylabel('Power (log(arb. unit))');
    %draw some vertical lines for better understanding of the plot
    fRange = [4,8,12,30,60,100];
    yline(0,':','LineWidth',2); %boundary for the x2/x1 (for actual oscillatory component
    xline(log(fRange),'--');
    text(log(fRange),ones(1,length(fRange))*3,cellstr(strtrim(num2str(fRange'))));
    xlim([0,ceil(log(max(fRange)))]);
    ylim([min(ylim),max(ylim)+1.1]);
    title(E_ft_laplace{i}.label{targetChannelIndex});
    pbaspect([1 1 1]);
    %legend({'Oscillatory'},'Location','SouthWest');
    hold off;
    drawnow;
    linkaxes(ax_all);
end
