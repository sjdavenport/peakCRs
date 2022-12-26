function out = latconvCRmean(data, top, use_stat_Lambda)
% latCR(data, top) finds confidence intervals for the location of
% maxima of the mean of data on a lattice
% -------------------------------------------------------------------------
% ARGUMENTS
% data    an object of class field.
% top     the number of peaks to find CRs at.
% use_stat_Lambda - optional boolean flag indicating whether to use a
%         stationary estimate of Lambda (i.e. estimated using all the data 
%         rather than the data at a point or not). Default is to use the
%         data at the point.
%--------------------------------------------------------------------------
% OUTPUT
%   out - structure containing the following fields:
%         max_locs - matrix of location of maxima of the mean of data
%         MFTD - the monte carlo distribution of the locations of the local maxima
%         peakderiv2 - matrix of second derivatives of the mean at the maxima
%         npeaks - number of peaks
%         peakvar - variance of the peak values
%         Lambda - cell array containing matrices of Lambda at the different peaks
%         Omega - cell array containing matrices of Omega at the different peaks
%         CR_lower - lower bounds of the confidence regions
%         CR_upper - upper bounds of the confidence regions
%         CR_volume - volume of the confidence regions
%--------------------------------------------------------------------------
% EXAMPLES
% % Compare to convCR
% lat_data = wfield([100,1], 100); 
% FWHM = 3; resadd = 21; params = ConvFieldParams(FWHM, resadd, 0);
% smooth_field = convfield(lat_data, params)
% out_latCR = latCRmean(smooth_field,1)
% peak_est_locs = {spacep_inv(lmindices(mean(smooth_field.field, 2), 1), resadd)};
% out_convCR = convCR(lat_data, FWHM, peak_est_locs)
% 
% % Stationary estimate of Lambda
% lat_data = wfield([100,1], 100); 
% FWHM = 3; resadd = 21; params = ConvFieldParams(FWHM, resadd, 0);
% smooth_field = convfield(lat_data, params)
% cutoff = ceil(4*FWHM2sigma(FWHM))*(resadd+1)
% cf = cut_field(smooth_field, cutoff);
% out = latconvCRmean(cf)
% out2 = latconvCRmean(cf, 1, 1)
% 
% increm = 0.001; smo = 5;
% peakspec = {[10,14], [21.5,26.5], [35,39]}; peakparams = {[1.5,2], [2,2], [2,1.5]};
% npeaks = length(peakspec);
% xvals = 1:increm:50; sigstore = peakgen( xvals, peakspec, peakparams, 1, smo);
% params = ConvFieldParams(FWHM, data_info.resadd);
% lat_data = sigstore' + statfield(Dim, nsubj, params, 2)*(1/sqrt(allvars(FWHM))); % Need to make noisegen take and output a Field!
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = data.D;
data_mean =  mean(data);
nsubj = data.fibersize;

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('top', 'var')
    top = 1;
end

if ~exist('use_stat_Lambda', 'var')
    use_stat_Lambda = 0;
end

%%  Main
%--------------------------------------------------------------------------
max_lat_locs = lmindices(data_mean.field, top, data.mask);
out.max_locs = data.xvals{1}(max_lat_locs); % returns a row vector!
[~,out.peakderiv2] = lat_derivs_dep(data_mean, out.max_locs);

out.npeaks = size(out.max_locs,2);
indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_derivs2_at_peaks = zeros(nsubj, D, out.npeaks); if D > 2; error('need to adjust this line'); end
indi_peak_vals = zeros(nsubj, out.npeaks);

if D == 1
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = data.field(max_lat_locs,subj);
        [indi_derivs_at_peaks(subj, :, :), indi_derivs2_at_peaks(subj, :, :)] = lat_derivs_dep(data(:,subj), out.max_locs);
    end
elseif D == 2
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, data(:,:,subj), Kernel, mask, truncation, xvals_vecs);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, data(:,:,subj), Kprime, mask, truncation, xvals_vecs);
    end
elseif D == 3
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, data(:,:,:,subj), Kernel, mask, truncation, xvals_vecs);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, data(:,:,:, subj), Kprime, mask, truncation, xvals_vecs);
    end
else
    error('Not coded yet')
end
out.peakvar = var(indi_peak_vals);

out.Lambda = cell(1, out.npeaks);

if use_stat_Lambda == 1
    lat_spacing = data.xvals{1}(2) - data.xvals{1}(1);
%     out.peakvar = m
    [~,~,Lambda] = est_smooth(data.field, data.mask, -1, lat_spacing); % Estimate Lambda from the data
    for I = 1:out.npeaks
        out.Lambda{I} = Lambda;
    end
else
    for I = 1:out.npeaks
        out.Lambda{I} = cov(indi_derivs_at_peaks(:, :, I));
        out.Omega{I} = cov(indi_derivs2_at_peaks(:, :, I));
        out.Delta{I} = (1/(nsubj-1))*sum((indi_derivs_at_peaks(:, :, I)-mean(indi_derivs_at_peaks(:,:,I))).*indi_derivs2_at_peaks(:, :, I));
        % use nsubj in the denominator as we haven't subtracted the mean
%         out.Delta{I} = (1/(nsubj))*sum((indi_derivs_at_peaks(:, :, I)-mean(indi_derivs_at_peaks)).*indi_derivs2_at_peaks(:, :, I));
        % (as one is mean zero!).
    end
end

out.cltSigmas = cell(1, out.npeaks);
out.muhessianinv = cell(1, out.npeaks);
out.dist95 = cell(1,out.npeaks);
out.asym95 = cell(1,out.npeaks);
alpha_quants = 0.05;
for I = 1:out.npeaks
    out.muhessianinv{I} = inv(reshape(out.peakderiv2(:,I), [D,D]));
    out.cltSigmas{I} = out.muhessianinv{I}*out.Lambda{I}*out.muhessianinv{I};
    if use_stat_Lambda == -1 %Gives you the option of estimating Lambda with or without
        Delta = 0;
    else
%         Delta = out.Delta{I};
        Delta = 0;
    end
    covmate = [out.Lambda{I}, Delta; Delta, out.Omega{I}];
    out.MFTD{I} = MFTD( out.peakderiv2(:,I), covmate, nsubj );
%     upper95dist = (1/2)*( - prctile(out.MFTD{I}, 100*alpha_quants/2));
    upper95asym = 1.96*sqrt(out.cltSigmas{I})/sqrt(nsubj);
    genpluspeak = out.max_locs(I) + out.MFTD{I}; 
    out.dist95{I} = [prctile(genpluspeak, 100*alpha_quants/2), prctile(genpluspeak, 100*(1-alpha_quants/2))];
    out.asym95{I} = [out.max_locs(I) - upper95asym, out.max_locs(I) + upper95asym];
end

end
