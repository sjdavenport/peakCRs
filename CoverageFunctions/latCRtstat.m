function [ out ] = latCRtstat(data, top)
% LATCRTSTAT( data, top ) creates confidence regions for the location of
% the maximum of Cohen's d.
%--------------------------------------------------------------------------
% ARGUMENTS
% data    an object of class field.
% top     the number of peaks to find CRs at.
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wfield([100,1], 50); 
% FWHM = 3; resadd = 21; params = ConvFieldParams(FWHM, resadd);
% smooth_field = convfield(lat_data, params)
% out = latCRtstat(smooth_field)
% out2 = tfieldCI(lat_data.field, lat_data.xvals, FWHM)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport Leave out out2
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = data.D;
[data_tstat, mu, sigma] =  mvtstat(data.field);
data_tstat = Field(data_tstat, 2); data_tstat.xvals = data.xvals;
nsubj = data.fibersize;

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('top', 'var')
    top = 1;
end

if ~exist('top', 'var')
    use_stat_Lambda = 0;
end

%%  Main
%--------------------------------------------------------------------------
max_lat_locs = lmindices(data_tstat.field, top, data.mask)';
out.max_locs = data.xvals{1}(max_lat_locs);
[~,tstatderiv2] = lat_derivs(data_tstat, out.max_locs);
out.CDderiv2 = num2cell(tstatderiv2/sqrt(nsubj));

out.npeaks = size(out.max_locs,2);
indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_derivs_at_peaks2 = zeros(nsubj, D, out.npeaks);
indi_peak_vals = zeros(nsubj, out.npeaks);

if D == 1
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = data.field(max_lat_locs,subj);
        indi_derivs_at_peaks(subj, :, :) = lat_derivs(data(:,subj), out.max_locs);
        indi_derivs_at_peaks2(subj, :, :) = lat_derivs(data(:,subj)./sigma, out.max_locs);
    end
% elseif D == 2
%     for subj = 1:nsubj
%         indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data(:,:,subj), Kernel,  mask, truncation, xvals_vecs);
%         indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data(:,:,subj), Kprime,  mask, truncation, xvals_vecs);
%     end
% elseif D == 3
%     for subj = 1:nsubj
%         indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data(:,:,:,subj), Kernel,  mask, truncation, xvals_vecs);
%         indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data(:,:,:, subj), Kprime, mask, truncation, xvals_vecs);
%     end
else
    error('Not coded yet')
end
out.peakvar = var(indi_peak_vals);

out.mu = mu(max_lat_locs);
outsigma = sigma(max_lat_locs);

out.sigma2 = outsigma.^2;

sigma2field = Field(sigma.^2, 2); sigma2field.xvals = data.xvals;
out.sigma2deriv = lat_derivs(sigma2field, out.max_locs); 

out.Lambda = cell(1, out.npeaks);
out.Lambdaprimeorig = cell(1, out.npeaks);
out.Lambdaprime = cell(1, out.npeaks);

out.Gamma = cell(1, out.npeaks);
out.derivvarorig = cell(1, out.npeaks);
out.derivvar = cell(1, out.npeaks);
out.cltSigmas = cell(1, out.npeaks);

for I = 1:out.npeaks
    demeaned_indi_peak_vals = indi_peak_vals - mean(indi_peak_vals,1);
    out.Gamma{I} =  mean(repmat(demeaned_indi_peak_vals(:, I),1,D).*indi_derivs_at_peaks(:, :, I),1)'; %Note only need to subtract the mean of one as this is covariance!!
    out.Lambda{I} = cov(indi_derivs_at_peaks(:, :, I)); 
    if I == 2
       a = 3
    end
    out.Lambdaprime{I} = cov(indi_derivs_at_peaks2(:, :, I)); 
    out.Lambdaprimeorig{I} = out.Lambda{I}/out.sigma2(I) - out.sigma2deriv(I)*out.Gamma{I}/out.sigma2(I)^2 + out.sigma2deriv(I)*out.sigma2deriv(I)'/(4*out.sigma2(I)^2);
    out.derivvarorig{I} = out.Lambdaprimeorig{I}*(1+(out.mu(I)^2/out.sigma2(I))*(nsubj/(nsubj-1))); %See Lemma 7.2 for the constant n/n-1 term. Is it in the right place or should it be outside of both brackets??
    out.derivvar{I} = out.Lambdaprime{I}*(1+(out.mu(I)^2/out.sigma2(I))*(nsubj/(nsubj-1)));
    out.cltSigmas{I} =  inv(out.CDderiv2{I})*out.derivvar{I}*inv(out.CDderiv2{I});
    upper95asym = 1.96*sqrt(out.cltSigmas{I})/sqrt(nsubj);
    out.asym95{I} = [out.max_locs(I) - upper95asym, out.max_locs(I) + upper95asym];
end 
 
end