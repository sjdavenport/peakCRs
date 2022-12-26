function [ ellval_store, maxlochist, coverage, bonf_coverage, Lambda_mate, peakderiv2mate,...
    cltSigmasmate] = convCIcoverage( data_info, niters, field_type, sim_type )
% convCIcoverage
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data_info   a data descriptor variable  
%  niters      the number of iterations
%  field_type  the type of field to simulate
%  sim_type    the type of simulation to run
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('sim_type', 'var')
    sim_type = 'N';
end

if ~exist('field_type', 'var')
    field_type = 'm';
end

% if strcmp(field_type, 'm')
%     meanonfinelat = data_info.meanonfinelat;
% end

meanfn = data_info.meanfn;
lat_std = data_info.lat_std;
FWHM = data_info.FWHM;
Dim = data_info.dim;
D = data_info.D;
npeaks = data_info.npeaks;

if ~isfield(data_info, 'xvals_vecs')
    data_info.xvals_vecs = {1:length(lat_signal)};
end
% xvalues_at_voxels = xvals2voxels(data_info.xvals_vecs);
nsubj = data_info.nsubj;

true_locs = data_info.maxlocs;

ellval_store = zeros(npeaks, niters);
maxlochist = zeros(data_info.D, npeaks, niters);

alpha_quants = [0.2, 0.1, 0.05];
nquants = length(alpha_quants);
alpha_quants = [alpha_quants, alpha_quants/npeaks];
coverage = zeros(npeaks,length(alpha_quants));
bonf_coverage = zeros(1, length(alpha_quants)/2);

% Load in the variance correction
global PIloc
load([PIloc,'Variance/storevars'], 'allvars')

Lambda_mate = zeros(D,D,npeaks,niters);
peakderiv2mate = zeros(D,D,npeaks,niters);
cltSigmasmate = zeros(D, D, npeaks,niters);

for L = 1:niters
    modul(L,1)
    if strcmp(sim_type, 'N') || strcmp(sim_type, 'Normal')
        % Scale by sqrt(allvars(FWHM)) to ensure the resulting fields are
        % variance 1.
        noise = lat_std*wfield(Dim, nsubj)*(1/sqrt(allvars(FWHM)));
    elseif strcmp(sim_type, 'tnoise')
        load([PIloc,'Variance/storevars_t3'], 'allvars')
        noise = lat_std*wfield(Dim, nsubj, 'T', 3)*(1/sqrt(allvars(FWHM)));
    elseif strcmp(sim_type, 'nonstat')
        rng(109)
        voxmap = randsample(dim(1), dim(1));
        noise = lat_std*wfield(Dim, nsubj, 'T', 3)*(1/sqrt(allvars(FWHM)));
    elseif strcmp(sim_type, 'chi2')
        noise = ((lat_std*wfield(Dim, nsubj)*(1/sqrt(allvars(FWHM)))).^2 -1)*(1/sqrt(2));
%     elseif strcmp(field_type, 'nonstat')
%         params = ConvFieldParams(FWHM, data_info.resadd, 0);
%         load([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatvar.mat'], 'var_est')
%         load([RFTboxloc, 'Random_Field_Generation/Nonstatnoisegeneration/nonstatvariables.mat'], 'Sigma', 'Sigmavars')
%         nvox = 10;
%         lat_data = zeros(nsubj,nvox);
%         for J = 1:nsubj
%             lat_data(J,:) = mvnrnd(zeros(1,nvox)', Sigma)./Sigmavars';
%         end
%         lat_data = lat_data';
%         lat_data = Field(lat_data, true(size(lat_data,1),1) );
%         lat_data = lat_signal' + lat_std*convfield(lat_data,params)./sqrt(var_est(FWHM-2,:))'; % Need to make noisegen take and output a Field!
    else
        error('This sim type has not been stored')
    end
    if isfield(data_info, 'xvals_lat')
        noise.xvals = data_info.xvals_lat;
    end
    
    peak_est_locs = cell(1, npeaks);
    for K = 1:npeaks
        peak_est_locs{K} = true_locs(:,K);
    end
    
    %     if strcmp(field_type, 'm')
    %         if D == 1
    %             peak_est_locs = npeaks;
    %         else
    %             peak_est_locs = true_locs;
    %         end
    %     else
    %         peak_est_locs = cell(1, npeaks);
    %         for K = 1:npeaks
    %             peak_est_locs{K} = true_locs(:,K);
    %         end
    %     end
    
    % Find the local maxima
    resadd_init = 1;
    if FWHM < 10 % This extra bit is only really needed for low FWHM
        if ~strcmp(field_type, 'm') && D > 1
            params = ConvFieldParams([FWHM, FWHM], resadd_init);
            f = convfield_t(noise + data_info.Sig, params);
            new_peak_est_locs = cell(1, length(peak_est_locs));
            %         for I = 1:length(peak_est_locs)
            for I = 1:length(peak_est_locs)
                round_peak_est_locs = (resadd_init+1)*round(peak_est_locs{I});
                addition = (resadd_init+1)*9;
                subset_f = f(round_peak_est_locs(1) - addition:round_peak_est_locs(1) + addition, ...
                    round_peak_est_locs(2) - addition:round_peak_est_locs(2) + addition);
                lms = lmindices(subset_f.field, 1);
                new_peak_est_locs{I} = xvaleval(lms, subset_f.xvals);
            end
            %         new_peak_est_locs = peak_est_locs;
            peak_est_locs = new_peak_est_locs;
        end
    end
    
    if strcmp(field_type, 'm')
%         out = convCR(noise, FWHM, meanfn, meanonfinelat.field, peak_est_locs, data_info.Sig);
        out = convCR(noise, FWHM, peak_est_locs, meanfn);
    else
        out = convCR_t(noise, FWHM, peak_est_locs, meanfn);
    end
    
    %     lat_data = data_info.lat_signal' + data_info.lat_std*smooth_noise; % Need to make noisegen take and output a Field!
    %     out2 = latCRtstat(smooth_data)
    %     out.cltSigmas{1}
    %     out2.cltSigmas{1}
    
    % Calculate the peak matching between the different peaks
    true_locs
    out.max_locs
    if npeaks > 1
        if D == 1
            obs2true_perm = match_peaks( true_locs, out.max_locs );
        else
            % Matching occurs simply by initializing at the true peak
            % location (not to be done that way in theory but fine in
            % practice as the peaks in our simulations will be strongly
            % identifiable).
            obs2true_perm = 1:npeaks;
        end
    else
        obs2true_perm = 1;
    end
    
    % Record the peak locations (after matching!)
    maxlochist(:,:,L) = out.max_locs(:,obs2true_perm);
    
    % For each (matched) peak calculate the L-value!
    bonf_store = zeros(1, length(alpha_quants)/2);
    for J = 1:npeaks
        ellval_store(J,L) = inellipse(true_locs(:,J), nsubj*inv(out.cltSigmas{obs2true_perm(J)}), out.max_locs(:,obs2true_perm(J)));
        Lambda_mate(:,:,J,L) = out.Lambda{obs2true_perm(J)};
        peakderiv2mate(:,:,J,L) = reshape(out.peakderiv2{obs2true_perm(J)}, [D,D]);
        cltSigmasmate(:,:,J,L) = out.cltSigmas{obs2true_perm(J)};
        
        if strcmp(field_type, 'm')
            temp_ells = zeros(1, length(out.MFTD{J}));
            for MFTD_iter = 1:length(out.MFTD{J})
                temp_ells(MFTD_iter) = inellipse(zeros(data_info.D,1), nsubj*inv(out.cltSigmas{obs2true_perm(J)}), out.MFTD{J}(MFTD_iter,:)');
            end
            
            for K = 1:length(alpha_quants)
                purse = prctile(temp_ells, 100*(1-alpha_quants(K)));
                if ellval_store(J,L) < purse
                    coverage(J,K) = coverage(J,K) + 1;
                    if K > nquants
                        bonf_store(K-nquants) = bonf_store(K-nquants) + 1;
                    end
                end
            end
            for alpha_q_iter = 1:length(alpha_quants)/2
                if bonf_store(alpha_q_iter) == npeaks
                    bonf_coverage(alpha_q_iter) = bonf_coverage(alpha_q_iter) + 1;
                end
            end
        end
        
    end
end

% Obtain the average
coverage = coverage/niters;
bonf_coverage = bonf_coverage/niters;

end

