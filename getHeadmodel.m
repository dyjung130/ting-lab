function getHeadmodel(dirName,fileName,eegSensorTemplate,plotModel)
fileToRead = fullfile(dirName,fileName);
mri = ft_read_mri(fileToRead,'dataformat','nifti');
%convert the coordinate system to 'mni' - for the Bolt et al., 2024 paper
mri = ft_convert_coordsys(mri,'mni','feedback','no');
%Do you want to change the anatomical labels for the axes [Y, n]? Y
%What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
%What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
%What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
%Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), s(scanner origin), n(ot a landmark)? n

%segment it: voxels of MRI are classified in three different tissue types:
%'brain', 'skull', and 'scalp'
cfg           = [];
cfg.output    = {'brain'};%,'skull','scalp'};
cfg.spmmethod = 'old';
segmentedmri  = ft_volumesegment(cfg, mri);
%check data
ft_checkdata(segmentedmri, 'feedback', 'yes');
%mesh it
cfg = [];
cfg.tissue = {'brain'};%,'skull','scalp'};
cfg.numvertices = [3000,1000,2000];
cfg.downsample = 2;
mesh = ft_prepare_mesh(cfg,segmentedmri);
%% mesh repair (using the iso2mesh functions in the fieldtrip folder)
% source: https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_mriproc.m#L97
for ii = 1:length(mesh)
    [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos,mesh(ii).tri,'dup');
    [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos,mesh(ii).tri,'isolated');
    [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos,mesh(ii).tri,'deep');
    [mesh(ii).pos, mesh(ii).tri] = meshcheckrepair(mesh(ii).pos,mesh(ii).tri,'meshfix');
end
%% save segmentedmri segmented mri
cfg = [];
cfg.method = 'dipoli';%'openmeeg';%singleshell for display only? it doesn't work for EEG source localize
vol = ft_prepare_headmodel(cfg, mesh);
vol = ft_convert_units(vol, 'mm');% I think the unit is already in mm but just in case it is not mm, do conversion

saveName = strsplit(fileName,'.');
saveName= saveName{1};
save(fullfile('mri',strcat('mriModel_',saveName)),'vol','mesh','cfg');

if plotModel == 1
    figure;
    ft_plot_sens(eegSensorTemplate);
    hold on;
    ft_plot_headmodel(vol);%results in mm
    title(saveName,'interpreter','none');
end
end