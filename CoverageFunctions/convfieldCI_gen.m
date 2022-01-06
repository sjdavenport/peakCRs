function out = convfieldCI_gen(lat_data, xvals_vecs, Kernel, peak_est_locs, mask)
% CONVFIELDCI_GEN(lat_data, xvals_vecs, Kernel, top, Kernelderiv, Kernelderiv2, alpha)
% finds confidence intervals for the location of maxima of a mean-convolution field
% -------------------------------------------------------------------------
% ARGUMENTS
% lat_data     a Dim by nsubj array giving the observed data
% xvalues_at_voxels    an nvoxel length vector giving the x values at the voxels
% Kernel    a function handle giving the kernel with which to smooth
% peak_est_locs  a D by npeaks matrix giving the initial estimates of the 
%               location of the peaks. If this is instead an integer: top  
%               then the top number of maxima are considered and initial 
%               locations are estimated from the underlying data. If this 
%               is not specified then it is set to 1, i.e. only considering
%               the maximum. If D = 1 and you wish to specify multiple
%               peaks rather than a number of peaks then you need to begin
%               your input as [NaN, peakestloc1, peakestloc2, ...].
% mask      a 0-1 array with the same dimensions as the data giving a mask
%--------------------------------------------------------------------------
% OUTPUT
% out.
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 5;
% noise = wfield([10,1], 100)
% out = convfieldCI_gen(noise.field, 1:10, FWHM)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 4
    peak_est_locs = 1;
end
if isnumeric(Kernel)
    FWHM = Kernel;
%     ss = @(x)applyconvfield_gen(x, mask, Kernelsquared, xvals_vecs);
%     h = 0.00001;
    Kernel = @(x) GkerMV(x,FWHM);
    Kprime = @(x) GkerMVderiv(x,FWHM);
    Kprime2 = @(x) GkerMVderiv2(x,FWHM);
%     Kernelsquared = @(x) GkerMV(x,FWHM).^2;
%     Kprime = @(x) [(Kernel(x+h*[1,0,0]') - Kernel(x))/h, (Kernel(x+h*[0,1,0]') - Kernel(x))/h, (Kernel(x+h*[0,0,1]') - Kernel(x))/h]';
%     Kprime2 = @(x) [(Kprime(x+h*[1,0,0]') - Kprime(x))/h, (Kprime(x+h*[0,1,0]') - Kprime(x))/h, (Kprime(x+h*[0,0,1]') - Kprime(x))/h];
end
truncation = 0;

Ldim = size(lat_data);
D = length(Ldim) - 1;
data_mean =  mean(lat_data, (D+1));
if D == 1
    data_mean = data_mean';
end
nsubj = size(lat_data, (D+1));

if nargin < 5
    mask = true(size(data_mean));
end

out.max_locs = findconvpeaks_orig(data_mean, FWHM, peak_est_locs, 'Z', mask, xvals_vecs);
out.peakderiv2 = applyconvfield(out.max_locs, data_mean, Kprime2, mask, truncation, xvals_vecs);

out.npeaks = size(out.max_locs,2);
indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_peak_vals = zeros(nsubj, out.npeaks);

if D == 1
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data(:,subj)', Kernel, mask, truncation, xvals_vecs);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data(:,subj)', Kprime, mask, truncation, xvals_vecs);
    end
elseif D == 2
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data(:,:,subj), Kernel, mask, truncation, xvals_vecs);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data(:,:,subj), Kprime, mask, truncation, xvals_vecs);
    end
elseif D == 3
    for subj = 1:nsubj
        indi_peak_vals(subj, :) = applyconvfield(out.max_locs, lat_data(:,:,:,subj), Kernel, mask, truncation, xvals_vecs);
        indi_derivs_at_peaks(subj, :, :) = applyconvfield(out.max_locs, lat_data(:,:,:, subj), Kprime, mask, truncation, xvals_vecs);
    end
end
out.peakvar = var(indi_peak_vals);

out.Lambda = cell(1, out.npeaks);
out.cltSigmas = cell(1, out.npeaks);
out.muhessianinv = cell(1, out.npeaks);
for I = 1:out.npeaks
    out.Lambda{I} = cov(indi_derivs_at_peaks(:, :, I)); 
    out.muhessianinv{I} = inv(reshape(out.peakderiv2(:,I), [D,D]));
    out.cltSigmas{I} = out.muhessianinv{I}*out.Lambda{I}*out.muhessianinv{I};
end

end
