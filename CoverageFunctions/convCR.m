function out = convCR(lat_data, FWHM, meanfn, peak_est_locs, use_stat_Lambda)
% CONVCR(lat_data, params, mean_function, top)
% finds confidence intervals for the location of maxima of the mean of a
% convolution field
% -------------------------------------------------------------------------
% ARGUMENTS
% lat_data     unsmoothed noise
% FWHM
% meanfn
% top
%--------------------------------------------------------------------------
% OUTPUT
% out.
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Obtain the number of dimensions and subjects
D = lat_data.D;
nsubj = lat_data.fibersize;

% Default the number of peaks to find to 1 (for now!)
if ~exist('top', 'var')
    top = 1;
end

if ~exist('use_stat_Lambda', 'var')
    use_stat_Lambda = -1;
end

% Set the Kernel truncation
truncation = ceil(7*FWHM2sigma(FWHM));
% truncation = 0;
Kprime = @(x) GkerMVderiv(x,FWHM);
Kprime2 = @(x) GkerMVderiv2(x,FWHM);

% Calculate derivative functions for the mean
[mean_derivfn, mean_deriv2fn] = MVderiv( meanfn, lat_data.D );

pre_smoothed_noise_mean = mean(lat_data);

% meannoisefield = @(x) applyconvfield(x, pre_smoothed_noise_mean, FWHM, lat_data.mask, truncation, lat_data.xvals);
meannoisederiv2 = @(x) applyconvfield(x, pre_smoothed_noise_mean.field, Kprime2, lat_data.mask, truncation, lat_data.xvals);

% if dochi2 == 0 
%     noisefn = @(x) applyconvfield(x, pre_smoothed_noise_mean.field, FWHM, lat_data.mask, truncation, lat_data.xvals);
% else
%     noisefn = @(x) applyconvfield(x, pre_smoothed_noise_mean.field, FWHM, lat_data.mask, truncation, lat_data.xvals).^2;
% end
% [~, meannoisederiv2] = MVderiv( noisefn, lat_data.D );

fneval = @(x) applyconvfield(x, pre_smoothed_noise_mean.field, FWHM, lat_data.mask, truncation, lat_data.xvals) + meanfn(x);
out.npeaks = length(peak_est_locs);
lowerbounds = cell(1, out.npeaks);
upperbounds = cell(1, out.npeaks);
for I = 1:out.npeaks
    lowerbounds{I} = peak_est_locs{I} - 6;
    upperbounds{I} = peak_est_locs{I} + 6;
end
peak_est_locs = cell2mat(peak_est_locs);
out.max_locs = findlms( fneval, peak_est_locs, lowerbounds, upperbounds );
% fcp = findconvpeaks(pre_smoothed_noise_mean, FWHM, peak_est_locs, 'Z', ceil(8*FWHM2sigma(FWHM)), 0, 1, meanfn, meanonlat)
% findlms( fneval, top, top - , upperbounds, algorithm )
% resadd = 51;
% params = ConvFieldParams([FWHM, FWHM], resadd);
% params.kernel.truncation = ceil(8*FWHM2sigma(FWHM));
% finenoise = convfield(pre_smoothed_noise_mean, params);
% add_FWHM = 1.5;
% params = ConvFieldParams([add_FWHM, add_FWHM], resadd);
% params.kernel.truncation = ceil(8*FWHM2sigma(add_FWHM));
% max_val = meanfn([49.5,49.5]');
% finemean = convfield(Sig, params);
% finefield = (finemean + finenoise)*(1/4.974043165253320);
% lms = lmindices(finefield.field, 3);
% out.max_locs
% maxlocs = xvaleval(lms, finefield.xvals)

out.npeaks = size(out.max_locs,2);

for I = 1:out.npeaks
    out.peakderiv2{I} = meannoisederiv2(out.max_locs(:,I)) + reshape(mean_deriv2fn(out.max_locs(:,I)), [D^2,1]);
end

indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_derivs2_at_peaks = zeros(nsubj, D*D, out.npeaks);
indi_peak_vals = zeros(nsubj, out.npeaks);

if D == 1
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs2_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,subj), Kprime2, lat_data.mask, truncation,  lat_data.xvals);
    end
elseif D == 2
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs2_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,subj), Kprime2, lat_data.mask, truncation,  lat_data.xvals);
    end
elseif D == 3
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,:,:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs2_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,:,subj), Kprime2, lat_data.mask, truncation,  lat_data.xvals);
    end
    
    % Unnecessary as really need properties of the noise here
    % I.e. basically adding a constant which is subtracted when doing
    % estimation anyways!
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = indi_peak_vals(subj, :) + meanfn(out.max_locs);
        indi_derivs_at_peaks(subj, :, :) = indi_derivs_at_peaks(subj, :, :)' + mean_derivfn(out.max_locs);
        indi_derivs2_at_peaks(subj, :, :) = indi_derivs2_at_peaks(subj, :, :) + reshape(mean_deriv2fn(out.max_locs), [D^2,1]);
    end
end

% Add the mean

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
        if D == 2
            out.Omega{I} = out.Omega{I}([1,2,4], [1,2,4]); % Obtain the unique rows!
        end
        out.Delta{I} = 0; %We're mainly going to use stationary fields for which Delta = 0!
%         out.Delta{I} = (1/(nsubj-1))*sum((indi_derivs_at_peaks(:, :, I)-mean(indi_derivs_at_peaks(:,:,I))).*indi_derivs2_at_peaks(:, :, I));        
% use nsubj in the denominator as we haven't subtracted the mean
    end
end

out.cltSigmas = cell(1, out.npeaks);
out.muhessianinv = cell(1, out.npeaks);
out.dist95 = cell(1,out.npeaks);
out.asym95 = cell(1,out.npeaks);
alpha_quants = 0.05;
for I = 1:out.npeaks
    out.muhessianinv{I} = inv(reshape(out.peakderiv2{I}, [D,D]));
    out.cltSigmas{I} = out.muhessianinv{I}*out.Lambda{I}*out.muhessianinv{I};
    if use_stat_Lambda == -1 %Gives you the option of estimating Lambda with or without
        Delta = zeros(D, D*(D+1)/2);
    else
        Delta = zeros(D, D*(D+1)/2);
        Delta = out.Delta{I};
    end
    covmate = [out.Lambda{I}, Delta; Delta', out.Omega{I}];
    out.MFTD{I} = MFTD( reshape(out.peakderiv2{I}, [D,D]), covmate, nsubj );
    upper95dist = (1/2)*( - prctile(out.MFTD{I}, 100*alpha_quants/2));
    upper95asym = 1.96*sqrt(out.cltSigmas{I})/sqrt(nsubj);
    genpluspeak = out.max_locs(I) + out.MFTD{I};
    out.dist95{I} = [prctile(genpluspeak, 100*alpha_quants/2), prctile(genpluspeak, 100*(1-alpha_quants/2))];
    out.asym95{I} = [out.max_locs(I) - upper95asym, out.max_locs(I) + upper95asym];
end

end
