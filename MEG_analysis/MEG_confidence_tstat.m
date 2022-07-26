global where_davenpor
prefft = load([where_davenpor,'Data/MEG/hcp_sails_aal_timeseries2.mat']);

%%
unique_ids = unique(prefft.subj_id);
count_vals = zeros(1,length(unique_ids));
for I = 1:length(unique_ids)
    count_vals(I) = sum(prefft.subj_id == unique_ids(I));
end
weirdone = unique_ids(count_vals < 3);

time_series = cell(1, length(unique_ids) - 1);
counter = 0;

ids2use = setdiff(unique_ids, weirdone);
for I = 1:length(ids2use)
    id = ids2use(I);
    indices = find(prefft.subj_id == id);
    time_series{I} = prefft.source_data{indices(1)};
end

%%
resadd = 111;
nsubj = 79;
power_spectrum = welch( time_series{1}, 240, 240, resadd );
ps_fields = Field(zeros(length(power_spectrum.mask), nsubj), power_spectrum.mask);
ps_fields.xvals{1} = power_spectrum.xvals{1};

for I = 1:nsubj
    I
    power_spectrum = welch( time_series{I}, 240, 240, resadd, 1 );
    ps_fields.field(:,I) = power_spectrum.field;
end

log_ps_fields = ps_fields;
log_ps_fields.field = log10(ps_fields.field);

%%

plot(log_ps_fields.xvals{1},mean(log_ps_fields.field, 2))

%%
plot(log_ps_fields.xvals{1}, mean(log_ps_fields.field,2)./std(log_ps_fields.field,0,2))

%% CRs for fields
outCRtstat = latCRtstat(log_ps_fields, 2)

%% CRs for fields minus the average
meanf = mean(log_ps_fields.field,2);
log_ps_fields_minus_min = log_ps_fields;
log_ps_fields_minus_min.field = log_ps_fields.field - min(meanf);
outCRtstat_minus_min = latCRtstat(log_ps_fields_minus_min, 2)

%% CRs for fields minus the average
meanf = mean(log_ps_fields.field,2);
endpoint = 60;
index = find(log_ps_fields.xvals{1} == endpoint);
average_log_power = mean(meanf(1:index));
log_ps_fields_minus_average = log_ps_fields;
log_ps_fields_minus_average.field = log_ps_fields.field - average_log_power;
outCRtstat_minus_average = latCRtstat(log_ps_fields_minus_average, 2)

%% Plot Cohen's d
set(0,'defaultAxesFontSize', 15); %This sets the default font size.

CD = mean(log_ps_fields.field,2)./std(log_ps_fields.field,0,2);
global PIloc
clf
for endpoint = {[5,15], [0,60]}
    startindex = find(log_ps_fields.xvals{1} == endpoint{1}(1));
    endindex = find(log_ps_fields.xvals{1} == endpoint{1}(2));
    h(1) = plot(log_ps_fields.xvals{1}(startindex:endindex),CD(startindex:endindex), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
    xlabel('Frequency (Hz)')
    ylabel('Cohen''s d')
    h(2) = xline( outCRtstat.asym95{1}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{1}(2), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(2), '--', 'LineWidth', 2);
    if endpoint{1}(2) == 60
        title('Cohen''s d of log power')
        legend(h(1:2), 'Cohen''s d', '95% Peak Confidence Interval', 'Location', 'SE')
    else
        title('Zoomed in Cohen''s d of log power')
    end
    export_fig([PIloc, 'EEG/MEG_confidence/Figures/tstat', num2str(endpoint{1}(2)), '.pdf'], '-transparent')
    pause
end
[1.96*sqrt(outCRtstat.cltSigmas{1})/sqrt(79), outCRtstat.max_locs(1)]
[1.96*sqrt(outCRtstat.cltSigmas{2})/sqrt(79), outCRtstat.max_locs(2)]



%% Plot Cohen's d when subtracting the average
% Compute the Bonferroni corrected uncertainty around each peak
asym = CRuncertainty( outCRtstat_minus_average, nsubj, 0.95 );

set(0,'defaultAxesFontSize', 15); %This sets the default font size.

CD = mean(log_ps_fields_minus_average.field,2)./std(log_ps_fields_minus_average.field,0,2);
global PIloc
clf
% ,[0,60]
for both = [0]
    for endpoint = {[0,3],[0,60]}
        clf
        startindex = find(log_ps_fields_minus_average.xvals{1} == endpoint{1}(1));
        endindex = find(log_ps_fields_minus_average.xvals{1} == endpoint{1}(2));
        h = plot(log_ps_fields_minus_average.xvals{1}(startindex:endindex),CD(startindex:endindex), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
        xlabel('Frequency (Hz)')
        ylabel('Cohen''s d')
        if endpoint{1}(2) == 3
            LW = 1.5;
        else
            LW = 2;
        end
        h(2) = xline( asym{1}(1), '--', 'LineWidth', LW);
        xline( asym{1}(2), '--', 'LineWidth', LW);
        xline( asym{2}(1), '--', 'LineWidth', LW);
        xline( asym{2}(2), '--', 'LineWidth', LW);
        
        if endpoint{1}(2) == 60
            title('Cohen''s d of log power')
%             legend(h(1:2), 'Cohen''s d', '95% Peak Confidence Interval', 'Location', 'SE')
        else
            ylim([0.6, 1.3])
            title('Zoomed in Cohen''s d of log power')
        end
        ylimits = get(gca,'Ylim');
%         for J = 1:2
%             point = outCRtstat_minus_average.max_locs(J);
%             if endpoint{1}(2) == 5
%                 width = asym{J}(2) - asym{J}(1);
%                 fill([point-width/2,point-width/2, point+width/2, point+width/2], ...
%                     [ylimits(1), ylimits(2), ylimits(2), ylimits(1)], ....
%                     [1 1 1]*0.75,'linestyle','none');
%             end
%             if both
%                 xline( point, '--', 'LineWidth', 1, 'color', 'black');
%             end
%         end
        ylim(ylimits)
        uistack(h,'top');
        
        export_fig([PIloc, 'EEG/MEG_confidence/Figures/tstat', num2str(endpoint{1}(2)), '_both_', num2str(both), '.pdf'], '-transparent')
    end
end


%% Plot Cohen's d for the minus min
set(0,'defaultAxesFontSize', 15); %This sets the default font size.

CD = mean(log_ps_fields_minus_min.field,2)./std(log_ps_fields_minus_min.field,0,2);
global PIloc
clf
for endpoint = {[5,15], [0,60]}
    startindex = find(log_ps_fields_minus_min.xvals{1} == endpoint{1}(1));
    endindex = find(log_ps_fields_minus_min.xvals{1} == endpoint{1}(2));
    h(1) = plot(log_ps_fields_minus_min.xvals{1}(startindex:endindex),CD(startindex:endindex), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
    xlabel('Frequency (Hz)')
    ylabel('Cohen''s d')
    h(2) = xline( outCRtstat.asym95{1}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{1}(2), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(2), '--', 'LineWidth', 2);
    if endpoint{1}(2) == 60
        title('Cohen''s d of log power')
        legend(h(1:2), 'Cohen''s d', '95% Peak Confidence Interval', 'Location', 'SE')
    else
        %         ylim([
        title('Zoomed in Cohen''s d of log power')
    end
    export_fig([PIloc, 'EEG/MEG_confidence/Figures/tstat', num2str(endpoint{1}(2)), '.pdf'], '-transparent')
end


%% Plot Cohen's d for the minus average
set(0,'defaultAxesFontSize', 15); %This sets the default font size.

CD = mean(log_ps_fields_minus_average.field,2)./std(log_ps_fields_minus_average.field,0,2);
global PIloc
clf
for endpoint = {[5,15], [0,60]}
    startindex = find(log_ps_fields_minus_average.xvals{1} == endpoint{1}(1));
    endindex = find(log_ps_fields_minus_average.xvals{1} == endpoint{1}(2));
    h(1) = plot(log_ps_fields_minus_average.xvals{1}(startindex:endindex),CD(startindex:endindex), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
    xlabel('Frequency (Hz)')
    ylabel('Cohen''s d')
    h(2) = xline( outCRtstat.asym95{1}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{1}(2), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(1), '--', 'LineWidth', 2);
    xline( outCRtstat.asym95{2}(2), '--', 'LineWidth', 2);
    if endpoint{1}(2) == 60
        title('Cohen''s d of log power')
        legend(h(1:2), 'Cohen''s d', '95% Peak Confidence Interval', 'Location', 'SE')
    else
        title('Zoomed in Cohen''s d of log power')
    end
    export_fig([PIloc, 'EEG/MEG_confidence/Figures/tstat', num2str(endpoint{1}(2)), '.pdf'], '-transparent')
end
