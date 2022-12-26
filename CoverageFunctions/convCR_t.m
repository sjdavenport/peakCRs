function [ out, tcf ] = convCR_t(lat_data, FWHM, meanfn, peak_est_locs )
% CONVCR_t(lat_data, FWHM, meanfn, peak_est_locs )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   lat_data    unsmoothed noise
%   FWHM        the smoothing applied to the data
%   meanfn      the mean function: for inference in practice this should be
%               @(x)0 but it may be useful to set this to be non-zero for
%               the purpose of simulations
%   peak_est_locs  a cell array containing the initial estimates of the peak locations
%--------------------------------------------------------------------------
% EXAMPLES
% nvox = 10; nsubj = 100; lat_data = wfield([nvox,1], nsubj)
% peak_centre = 5; peakspec = {[3,7]}; peakparams = {[2,2]};
% [~, meanfn] = peakgen1D( 1:0.1:nvox, peakspec, peakparams, 1, 0.01);
% FWHM = 2; out = convCR_t(lat_data, FWHM, meanfn, {5});
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

xvals_vecs = lat_data.xvals;
mask = lat_data.mask;
D = lat_data.D;
nsubj = lat_data.fibersize;

% Default is to take the mean to be 0
if ~exist('meanfn', 'var')
    meanfn = @(x) 0;
end

% Default is to target the maximum
if ~exist('peak_est_locs', 'var')
    peak_est_locs = 1;
end

Kprime = @(x) GkerMVderiv(x,FWHM);
truncation = 0;

tcf = @(tval) applyconvfield_t( tval, lat_data.field, FWHM, mask, truncation, xvals_vecs, meanfn );
h = 0.0001;
if D == 1
    tcf_deriv = @(tval) (tcf(tval+h) - tcf(tval))/h;
    tcf_deriv2 = @(tval) (tcf_deriv(tval+h) - tcf_deriv(tval))/h;
elseif D == 2
    tcf_deriv = @(x) [(tcf(x+h*[1,0]') - tcf(x))/h, (tcf(x+h*[0,1]') - tcf(x))/h]';
    tcf_deriv2 = @(x) [(tcf_deriv(x+h*[1,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,1]') - tcf_deriv(x))/h];
elseif D == 3
    tcf_deriv = @(x) [(tcf(x+h*[1,0,0]') - tcf(x))/h, (tcf(x+h*[0,1,0]') - tcf(x))/h, (tcf(x+h*[0,0,1]') - tcf(x))/h]';
    tcf_deriv2 = @(x) [(tcf_deriv(x+h*[1,0,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,1,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,0,1]') - tcf_deriv(x))/h];
end

out.npeaks = length(peak_est_locs);
lowerbounds = cell(1, out.npeaks);
upperbounds = cell(1, out.npeaks);
for I = 1:out.npeaks
    lowerbounds{I} = peak_est_locs{I} - 4;
    upperbounds{I} = peak_est_locs{I} + 4;
end
peak_est_locs = cell2mat(peak_est_locs);
out.max_locs = findlms(tcf, peak_est_locs, lowerbounds, upperbounds);

indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_peak_vals = zeros(nsubj, out.npeaks);

if D == 1
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
    end
elseif D == 2
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
    end
elseif D == 3
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data.field(:,:,:,subj), FWHM, lat_data.mask, truncation,  lat_data.xvals);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data.field(:,:,:,subj), Kprime, lat_data.mask, truncation,  lat_data.xvals);
    end
    
    % Unnecessary as really need properties of the noise here
    % I.e. basically adding a constant which is subtracted when doing
    % estimation anyways!
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = indi_peak_vals(subj, :) + meanfn(out.max_locs);
        indi_derivs_at_peaks(subj, :, :) = indi_derivs_at_peaks(subj, :, :)' + mean_derivfn(out.max_locs);
    end
end
out.peakvar = var(indi_peak_vals);

[ ~, out.mu, outsigma ] = applyconvfield_t( out.max_locs, lat_data.field, FWHM, mask, truncation, xvals_vecs, meanfn );

out.sigma2 = outsigma.^2;
sigma2plush = zeros(out.npeaks,D); 
out.sigma2deriv = zeros(out.npeaks, D);
sb_vectors = eye(D);
for d = 1:D
    [ ~, ~, sigmaplushd ] = applyconvfield_t( out.max_locs+h*sb_vectors(:,d), lat_data.field, FWHM, mask, truncation, xvals_vecs, meanfn );
    sigma2plush(:, d) = sigmaplushd.^2;
end
out.sigma2deriv = (sigma2plush - out.sigma2)/h; %Note better to look at the derivative of sigma^2 
% rather than sigma as sigmahat^2 is an unbiased estimator of sigma^2! whereas the same is not true of sigmahat

% deriv_clt_var = out.Lambda./out.sigma.^2 - out.Gamma.*out.sigma2deriv./out.sigma.^4 + out.sigma2deriv.^2./out.sigma.^6;
% deriv_clt_var = deriv_clt_var.*(1 + out.mu.^2./out.sigma.^2);
% out.clt_std = abs(sqrt(deriv_clt_var)/out.CDderiv2);

out.Lambda = cell(1, out.npeaks);
out.Gamma = cell(1, out.npeaks);
% out.CDderiv2 = cell(1, out.npeaks);
out.derivvar = cell(1, out.npeaks);
out.cltSigmas = cell(1, out.npeaks);

for I = 1:out.npeaks
    demeaned_indi_peak_vals = indi_peak_vals - mean(indi_peak_vals,1);
    out.Gamma{I} =  mean(repmat(demeaned_indi_peak_vals(:, I),1,D).*indi_derivs_at_peaks(:, :, I),1)'; %Note only need to subtract the mean of one as this is covariance!!
    out.Lambda{I} = cov(indi_derivs_at_peaks(:, :, I));
    out.derivvar{I} = out.Lambda{I}/out.sigma2(I) - out.sigma2deriv(I,:)*out.Gamma{I}/out.sigma2(I)^2 + out.sigma2deriv(I,:)*out.sigma2deriv(I,:)'/(4*out.sigma2(I)^2);
    out.derivvar{I} = out.derivvar{I}*(1+(out.mu(I)^2/out.sigma2(I))*(nsubj/(nsubj-1))); %See Lemma 7.2 for the constant n/n-1 term. Is it in the right place or should it be outside of both brackets??
    %could rep out.mu(I)^2/out.sigma2(I) with the t-stat scaled of course
    %to see if it makes a difference???
    out.peakderiv2{I} = tcf_deriv2(out.max_locs(:,I))/sqrt(nsubj);
    out.cltSigmas{I} =  inv(out.peakderiv2{I})*out.derivvar{I}*inv(out.peakderiv2{I});
end
 
% [ ~, muminush, sigmaminush ] = tcfield( out.peak_locs-h, lat_data, xvalues_at_voxels, Kernel );
% out.muderiv2 = (muplush - 2*out.mu + muminush)/h^2;
% out.sigma2deriv2 = (sigmaplush - 2*out.sigma.^2 + sigmaminush)/h^2;
end

% 
% if isequal(size(peak_est_locs), [1,1])
%     top = peak_est_locs;
%     xvalues_at_voxels = xvals2voxels( xvals_vecs );
%     teval_lat = tcf(xvalues_at_voxels);
%     max_indices = lmindices(teval_lat, top, mask)'; %Note the transpose here! It's necessary for the input to other functions.
%     if D == 1
%         max_indices = max_indices';
%     end
%     top = length(max_indices);  
%     peak_est_locs = zeros(D, top);
%     for I = 1:D
%         peak_est_locs(I, :) = xvals_vecs{I}(max_indices(I,:));
%     end
% end
% 
% %Could use normal convolution methods above as just need the ones on the lattice
% % out.peak_locs = approx_peak_locs(lat_data, Kernel, xvalues_at_voxels, peak_est_locs);
% if D == 1 && isnan(peak_est_locs(1)) %This allows for multiple 1D peaks!
%     peak_est_locs = peak_est_locs(2:end);
% end
% npeaks = size(peak_est_locs, 2);
% 
% out.max_locs = zeros(D, npeaks);
% for peakI = 1:npeaks
% %     applyconvfield_gen(peak_est_locs(:, peakI), lat_data, Kprime, xvals_vecs )
% %     field_deriv(peak_est_locs(:, peakI))
%     out.max_locs(:, peakI) = NewtonRaphson(tcf_deriv, peak_est_locs(:, peakI), tcf_deriv2);
% end
