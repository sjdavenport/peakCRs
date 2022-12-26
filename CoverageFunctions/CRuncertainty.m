function [ asym, montecarlo ] = CRuncertainty( out, nsubj, quantile, do_bonf )
% CRuncertainty( out, nsubj, quantile, do_bonf )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  out    the output from convCR or convCR_t
%  nsubj  the number of subjects
%  quantile   the desired quantile
% Optional
%  do_bonf    0/1 whether to apply a bonferroni correction
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'do_bonf', 'var' )
   % Default value
   do_bonf = 1;
end

if ~exist( 'quantile', 'var' )
   % Default value
   quantile = 0.95;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
%% Compute the confidence intervals
npeaks = length(out.max_locs);
asym = cell(1,npeaks);

if isfield(out, 'dist95')
    montecarlo = cell(1,npeaks);
else
    montecarlo = NaN;
end
alpha_quants = (1-quantile)/npeaks;

% If no Bonferroni correction use 1 peak
if ~do_bonf
    npeaks = 1;
end

% Compute the uncertainty
for I = 1:2
    upper95asym = abs(norminv((1-quantile)/2/npeaks))*sqrt(out.cltSigmas{I})/sqrt(nsubj);
    if isfield(out, 'dist95')
        genpluspeak = out.max_locs(I) + out.MFTD{I};
        montecarlo{I} = [prctile(genpluspeak, 100*alpha_quants/2), prctile(genpluspeak, 100*(1-alpha_quants/2))];
    end
    asym{I} = [out.max_locs(I) - upper95asym, out.max_locs(I) + upper95asym];
    
    disp('Peak Location:')
    out.max_locs(I)
    disp('Asymptotic Uncertainty:')
    (asym{I}(2) - asym{I}(1))/2
    disp('Asymptotic Interval:')
    fprintf('(%.3f,%.3f)\n',  asym{I}(1), asym{I}(2))
    if isfield(out, 'dist95')
        disp('Monte Carlo Uncertainty:')
        (montecarlo{I}(2) - montecarlo{I}(1))/2
        disp('Monte Carlo Interval:')
        fprintf('(%.3f,%.3f)\n',  montecarlo{I}(1), montecarlo{I}(2))
    end
end

end

