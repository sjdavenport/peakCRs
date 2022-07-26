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

% log_ps_fields.field = log_ps_fields.field - 24.85;
%%
meanf = mean(log_ps_fields.field, 2);
endpoint = 60;
index = find(log_ps_fields.xvals{1} == endpoint);
average_log_power = mean(meanf(1:index));

%%
plot(log_ps_fields.xvals{1},mean(log_ps_fields.field, 2))

%%
plot(mean(log_ps_fields.field)./std(log_ps_fields.field))

%%
outCRmean = latCRmean(log_ps_fields, 2)

%% Plot the noise!
meanf = mean(log_ps_fields.field, 2);
nvox = 100;
for I = 1:6
    noise = log_ps_fields.field(:,I) - meanf;
    subplot(2,1,1); plot(meanf(1:nvox))
    subplot(2,1,2); plot( noise(1:nvox));
    pause
end

%%
dist = CRuncertainty( outCRmean, nsubj, 0.95 );

%% Plot the fields!
set(0,'defaultAxesFontSize', 15); %This sets the default font size.

meanf = mean(log_ps_fields.field, 2);
global PIloc
for both = [0]
    for endpoint = [3,60]
        clf
        index = find(log_ps_fields.xvals{1} == endpoint);
        h(1) = plot(log_ps_fields.xvals{1}(1:index),meanf(1:index), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
        for I = 1:6
            %     clf
            field = log_ps_fields.field(:,I+2);
            h(I+1) = plot( log_ps_fields.xvals{1}(1:index), field(1:index), 'color', def_col('blue'), 'Linewidth', 1);
            %     pause
        end
        ylimits = get(gca,'Ylim');
        for J = 1:2
            point = outCRmean.max_locs(J);
%             if endpoint == 3
%                 width = dist{J}(2) - dist{J}(1);
%                 fill([point-width/2,point-width/2, point+width/2, point+width/2], ...
%                     [ylimits(1), ylimits(2), ylimits(2), ylimits(1)], ....
%                     [1 1 1]*0,'linestyle','none')
% %                 fill([point-width/2,point-width/2, point+width/2, point+width/2], ...
% %                     [ylimits(1), ylimits(2), ylimits(2), ylimits(1)], ....
% %                     [1 1 1]*0.75,'linestyle','none')
%             end
%             if both
%                 xline( point, '--', 'LineWidth', 1, 'color', 'black')
%             end
        end
        ylim(ylimits)
        uistack(h,'top');
        
        h(I+1) = xline( dist{1}(1), '--', 'LineWidth', 1.5);
        xline( dist{1}(2), '--', 'LineWidth', 1.5);
        xline( dist{2}(1), '--', 'LineWidth', 1.5);
        xline( dist{2}(2), '--', 'LineWidth', 1.5);
        xlabel('Frequency (Hz)')
        ylabel('log( Power )')
        if endpoint == 5
            yticks([26, 26.5]);
            title('Zoomed in average log power spectrum')
            %         xline(2.26, 'LineWidth', 2)
            %         xline(2.32, 'LineWidth', 2)
        else
            title('Average log power spectrum')
            legend(h(1:2), 'Mean over Subjects', 'Individual Subjects')
        end
        export_fig([PIloc, 'EEG/MEG_confidence/Figures/mean', num2str(endpoint), '_both_', num2str(both), '.pdf'], '-transparent')
    end
end
%% Plot the fields minus the average!
set(0,'defaultAxesFontSize', 15); %This sets the default font size.

log_ps_fields_minus_average = log_ps_fields;
log_ps_fields_minus_average.field = log_ps_fields.field - average_log_power;

meanf = mean(log_ps_fields_minus_average.field, 2);
global PIloc
for endpoint = [10,60]
    clf
    index = find(log_ps_fields_minus_average.xvals{1} == endpoint);
    h(1) = plot(log_ps_fields_minus_average.xvals{1}(1:index),meanf(1:index), 'color', def_col('blue'), 'Linewidth', 2.5); hold on;
    for I = 1:6
        %     clf
        field = log_ps_fields_minus_average.field(:,I+2);
        h(I+1) = plot( log_ps_fields_minus_average.xvals{1}(1:index), field(1:index), 'color', def_col('blue'), 'Linewidth', 1);
        %     pause
    end
    legend(h(1:2), 'Mean over Subjects', 'Individual Subjects')
    xlabel('Frequency (Hz)')
    ylabel('log( Power ) - log( Power_{av} )')
    if endpoint == 10
        %         yticks([26, 26.5]);
        title('Zoomed in average log power spectrum')
        %         xline(2.26, 'LineWidth', 2)
        %         xline(2.32, 'LineWidth', 2)
    else
        title('Average log power spectrum')
    end
    export_fig([PIloc, 'EEG/MEG_confidence/Figures/mean', num2str(endpoint), '.pdf'], '-transparent')
end

%% DEP

%%
plot(mean(ps_fields.field, 2))

%%
plot(mean(ps_fields.field)./std(ps_fields.field))

%% Plot the noise!
meanf = mean(ps_fields.field, 2);
nvox = 100;
for I = 1:79
    noise = ps_fields.field(:,I) - meanf;
    subplot(2,1,1); plot(meanf(1:nvox))
    subplot(2,1,2); plot( noise(1:nvox));
    pause
end

%%
plot(10*log10(ps_fields.field(:,I)))

%%
plot(mean(ps_fields.field, 2))
%% Wierdly here the dist95 is less than the asym95?? How odd!!
outCRmean = latCRmean(ps_fields./10^26)

%%
outCRtstat = latCRtstat(ps_fields./10^26)

