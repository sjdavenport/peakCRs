function out = latCRmeanND(data, top, use_stat_Lambda)
% latCRmeanND(data, top) finds confidence intervals for the location of
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
% % Compare to convCR to ensure they have the same coverage
% lat_data = wfield([20,20], 50); 
% FWHM = 6; resadd = 7; params = ConvFieldParams([FWHM,FWHM], resadd, 0);
% smooth_field = convfield(lat_data, params)
% out_latCR = latCRmeanND(smooth_field,1);
% peak_est_locs = {spacep_inv(lmindices(mean(smooth_field.field, 3), 1), resadd)};
% out_convCR = convCR(lat_data, FWHM, peak_est_locs);
% out_latCR.cltSigmas{1}
% out_convCR.cltSigmas{1}
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
out.max_locs = zeros(size(max_lat_locs));
for d = 1:D
    out.max_locs(d,:) = data.xvals{d}(max_lat_locs(d,:));
end
[~,out.peakderiv2] = lat_derivs(data_mean, out.max_locs);

% Obtain the number of peaks
out.npeaks = size(out.max_locs,2);

% Initalize arrays to store the first and second derivatives at the
% different peaks
indi_derivs_at_peaks = zeros(nsubj, D, out.npeaks);
indi_derivs2_at_peaks = zeros(nsubj, D, D, out.npeaks); 
indi_peak_vals = zeros(nsubj, out.npeaks);

% Convert the maximum locations to cells for multiple indexing
maxlocsascell = cell(1,D);
for d = 1:D
    maxlocsascell{d} = max_lat_locs(d, :);
end
    
for subj = 1:nsubj
    indi_peak_vals(subj, :) = diag(data.field(maxlocsascell{:},subj));
    var_index = repmat({':'}, 1, D); var_index{D+1} = subj;
    [indi_derivs_at_peaks(subj, :, :), indi_derivs2_at_peaks(subj, :, :, :)] = lat_derivs(data(var_index{:}), out.max_locs);
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
    vech_indexing = triu(vech(1:(D*(D+1)/2)));
    for I = 1:out.npeaks
        out.Lambda{I} = cov(indi_derivs_at_peaks(:, :, I));
        out.Omega{I} = zeros(D*(D+1)/2, D*(D+1)/2);
        for J = 1:(D*(D+1)/2)
            Jindex = convind(find(vech_indexing == J), [D,D]);
            if D == 1
                Jindex = [1,1];
            end
            for K = 1:(D*(D+1)/2)
                Kindex = convind(find(vech_indexing == K), [D,D]);
                if D == 1
                    Kindex = [1,1];
                end
                out.Omega{I}(J,K) = ...
                    (1/(nsubj-1))*sum((indi_derivs2_at_peaks(:, Jindex(1), Jindex(2), I)-...
                    mean(indi_derivs2_at_peaks(:,Jindex(1),Jindex(2),I), 1)).*indi_derivs2_at_peaks(:, Kindex(1), Kindex(2), I));
            end
        end
        
        % For now only compute Delta = cov(\nabla Y, V(\nabla^2 Y)) if D = 1 
        if D > 1
            out.Delta{I} = zeros(D, D*(D+1)/2);
        else
            out.Delta{I} = (1/(nsubj-1))*sum((indi_derivs_at_peaks(:, :, I)-mean(indi_derivs_at_peaks(:,:,I))).*indi_derivs2_at_peaks(:, :, I));
        end
    end
end

out.cltSigmas = cell(1, out.npeaks);
out.muhessianinv = cell(1, out.npeaks);
out.dist95 = cell(1,out.npeaks);
out.asym95 = cell(1,out.npeaks);
alpha_quants = 0.05;

for I = 1:out.npeaks
    out.muhessianinv{I} = inv(out.peakderiv2(:,:,I));
    out.cltSigmas{I} = out.muhessianinv{I}*out.Lambda{I}*out.muhessianinv{I};
    if use_stat_Lambda == -1 %Gives you the option of estimating Lambda with or without
        Delta = zeros(D, D*(D+1)/2);
        Delta = 0;
    else
        Delta = zeros(D, D*(D+1)/2);
        Delta = out.Delta{I};
    end
    covmate = [out.Lambda{I}, Delta; Delta', out.Omega{I}];
    out.MFTD{I} = MFTD( out.peakderiv2(:,:,I), covmate, nsubj );    
    upper95asym = 1.96*sqrt(out.cltSigmas{I})/sqrt(nsubj);
    genpluspeak = out.max_locs(I) + out.MFTD{I}; 
    out.dist95{I} = [prctile(genpluspeak, 100*alpha_quants/2), prctile(genpluspeak, 100*(1-alpha_quants/2))];
    out.asym95{I} = [out.max_locs(I) - upper95asym, out.max_locs(I) + upper95asym];
end

end
