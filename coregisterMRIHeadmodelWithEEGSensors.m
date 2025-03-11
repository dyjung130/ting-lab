function coregisterMRIHeadmodelWithEEGSensors(dirName,fileName,eegSensorTemplate,plotModel)
%for subject i = 9
fileToRead = fullfile(dirName,fileName);
mri = ft_read_mri(fileToRead,'dataformat','nifti');
%convert the coordinate system to 'mni' - for the Bolt et al., 2024 paper
mri = ft_convert_coordsys(mri,'mni','feedback','no');
%Do you want to change the anatomical labels for the axes [Y, n]? Y
%What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
%What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
%What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
%Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), s(scanner origin), n(ot a landmark)? n

%segment it!
cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri);
%save segmentedmri segmented mri
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);
vol = ft_convert_units(vol, 'mm');% I think the unit is already in mm but just in case it is not mm, do conversion

saveName = strsplit(fileName,'.');
saveName= saveName{1};
save(fullfile('mri',strcat('mriModel_',saveName)),'vol','cfg');

if plotModel == 1
    figure;
    ft_plot_sens(eegSensorTemplate);
    hold on;
    ft_plot_headmodel(vol);%results in mm
    title(saveName,'interpreter','none');
end
end