function viewcr( crout, smooth_mean, crtype, quant )
% VIEWCR( crout, mask, crtype, slice ) plots the confidence regions for the
% peaks of a random field.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   crout    - a struct containing the output of the function convCR.
%   smooth_mean - a struct containing the output of the function convfield.
% Optional
%   crtype - a string specifying the type of confidence region to plot 
%            (either 'asymptotic' or 'MC'). If not provided, the default 
%              value is 'MC'.
%   quant - a scalar in the range [0, 1] specifying the confidence level 
%               of the confidence region. If not provided, the default value is 0.95.
%--------------------------------------------------------------------------
% OUTPUT
% No output
%--------------------------------------------------------------------------
% EXAMPLES
% nsubj = 50; dim = [15,15];
% lat_data = wfield(dim, nsubj);
% unsmooth_sig = peakgen(0.25, 1, 4, dim);
% lat_data.field = lat_data.field + unsmooth_sig;
% FWHM = 2; out = convCR(lat_data, FWHM, {(dim/2)'})
% resadd = 3; params = ConvFieldParams( [FWHM, FWHM], resadd );
% smoothfield = convfield(lat_data, params);
% viewcr( out, mean(smoothfield) )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
mask = smooth_mean.mask;
color = zeros([size(mask), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-mask);

if ~exist('crtype', 'var')
    crtype = 'MC';
end

if ~exist('quant', 'var')
    quant = 0.95;
end

if strcmp(crtype, 'asym') || strcmp(crtype, 'asymptotic')
    useMFTD = 0;
else
    if strcmp(crtype, 'MC') || strcmp(crtype, 'MCMC ') || strcmp(crtype, 'MonteCarlo')
        useMFTD = 1;
    else
        error('Only Monte carlo and asymptotic confidence intervals are available')
    end
end

% if useMFTDcov
%     covmateinv = inv(cov(crout.MFTD{I}));
% else
%     covmateinv = nsubj*inv(crout.cltSigmas{I});
% end

%%  Main Function Loop
%--------------------------------------------------------------------------
npeaks = length(crout.MFTD);
for I = 1:npeaks
    covmateinv = smooth_mean.fibersize*inv(crout.cltSigmas{I});

    temp_ells = zeros(1, length(crout.MFTD{I}));
    
    % Calculate CR quantile
    if useMFTD
        for MFTD_iter = 1:length(crout.MFTD{I})
            temp_ells(MFTD_iter) = inellipse(zeros(2,1), covmateinv, crout.MFTD{I}(MFTD_iter,:)');
        end
        chi2quant = prctile(temp_ells, 100*quant);
    else
        chi2quant = chi2inv(quant, 2);
    end
    
    crdata{I} = zeros(smooth_mean.masksize);
    for K = 1:length(smooth_mean.xvals{1})
        K
        for J = 1:length(smooth_mean.xvals{2})
            point = [smooth_mean.xvals{1}(K), smooth_mean.xvals{2}(J)]';
            [~, crdata{I}(K, J)] = inellipse(point, covmateinv, crout.max_locs(:,I), chi2quant);
        end
    end
    colored_slice{I} = zeros([size(fliplr(crdata{I})'), 3]);
    colored_slice{I}(:,:,1) = fliplr(crdata{I})';
end

subplot(1,2,1)
surf(smooth_mean.field)

subplot(1,2,2)
imagesc(fliplr(smooth_mean.field)')
hold on
for I = 1:npeaks
    im2 = imagesc(colored_slice{I});
    set(im2,'AlphaData',fliplr(crdata{I})');
    hold on
end
color = zeros([size(fliplr(smooth_mean.mask)'), 3]);
im2 = imagesc(color);
set(im2,'AlphaData',1-fliplr(smooth_mean.mask)');
axis off

end

