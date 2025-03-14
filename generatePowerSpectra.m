%input,
%datEEG: EEG data in FieldTrip format
%datLatency: Time stamps of the data to be analyzed; e.g., [1 100] in seconds
%freqRange: Frequency range; e.g., [1,100] in Hertz
function generatePowerSpectra(datEEG,datLatency,freqRange)
parpool("processes",4);%setup parallel here.
%% FOOOF for all channels
%normalize data first
datEEGz = [];% Z-scored EEG data

for i = 1:length(datEEG)
    datEEGz{i} = datEEG{i};
    datEEGz{i}.trial = cellfun(@(x) zscore(x,[],'all'),datEEGz{i}.trial,'uni',0);
end

set(0,'DefaultFigureWindowStyle','docked')

figure;
ax_all = [];
parfor targChanInd = 1:length(datEEG{1}.label)

    datOsc = [];%oscillatory power

    for i = 1:length(datEEG)
        %shorten data first (for testing)
        cfg = [];
        cfg.latency = datLatency;
        cfg.channel = targChanInd;%channe of interest
        datEEGzs = ft_selectdata(cfg,datEEGz{i});%segmented EEGz data

        cfg = [];
        cfg.method    = 'mtmfft';
        cfg.output    = 'fooof_aperiodic';
        cfg.pad='nextpow2';
        cfg.tapsmofrq = 4;
        cfg.foilim    = freqRange; %Hz
        E_1overf = ft_freqanalysis(cfg, datEEGzs);%aperiodic component (fractal)

        cfg.output = 'pow';
        E_orig = ft_freqanalysis(cfg,datEEGzs);%original spectrum

        %subtract the 1/f component from the power spectrum
        cfg =  [];
        cfg.parameter = 'powspctrm';
        cfg.operation = 'x2./x1';
        datOsc{i} = ft_math(cfg,E_1overf,E_orig);%oscillations
    end

    powSpectra = cellfun(@(x) x.powspctrm, datOsc, 'uni',0);
    powSpectra = log(vertcat(powSpectra{:}));%concatenate (subjects x freqpow)
    %
    meanSpectrum = mean(powSpectra,1);
    stdSpectrum = std(powSpectra,[],1);
    semSpectrum = stdSpectrum/sqrt(size(powSpectra,1)-1);

    ax_all(targChanInd) = subplot(4, ceil(length(datEEG{i}.label)/4),targChanInd);
    hold on;
    %plot(log(E_osc.freq),log(E_osc.powspctrm),'LineWidth',2)
    shadedErrorBar(log(datOsc{1}.freq),(meanSpectrum),(semSpectrum),'lineprop','r','patchSaturation',0.12)
    %xlabel('Frequency (log(Hz))');
    %ylabel('Power (log(arb. unit))');
    %draw some vertical lines for better understanding of the plot
    fRange = [4,8,12,30,60,100];
    yline(0,':','LineWidth',2); %boundary for the x2/x1 (for actual oscillatory component
    xline(log(fRange),'--');
    %text(log(fRange),ones(1,length(fRange)),cellstr(strtrim(num2str(fRange'))),'FontSize',14);
    xlim([log(min(fRange)),ceil(log(max(fRange)))]);
    ylim([min(ylim),1.1]);
    title(datEEG{i}.label{targChanInd});
    pbaspect([1 1 1]);
    set(gca,'fontsize',16);
    set(gca,'xtick',log(fRange),'xticklabel',num2str(fRange'),'XTickLabelRotation',45)
    hold off;
    drawnow;

end
linkaxes(ax_all);