function data_info = init_data_info1D( FWHM, nsubj, x, peakspec, peakparams, peakheights, smo)
% INIT_DATA_INFO1D initializes the data_info framework for running coverage
% testing in 1D
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM    the smoothness of the noise
%  nsubj   the number of subjects
%  x       the vector at which to evaluate the signal
%  peakspec   a cell array each entry being a 2D vector with start and
%             ending vertices for the peak to occur
%  peakparams  a cell array each entry being a 2D vector specifying the
%              peak parameters (used to generate the peak from the beta pdf)
%  peakheights  a vector giving the height of each peak
% Optional
%  smo   the amount of additional smoothing to apply to the data  (not
%        implemented yet!)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
Sig = peakgen1D( x, peakspec, peakparams, peakheights, smo ); add_FWHM = 1.5;
resadd = 0; params = ConvFieldParams([add_FWHM,add_FWHM], resadd);
meanonlat = convfield(Sig, params);
params.resadd = 1; meanonfinelat = convfield(Sig, params);
truncation = 4*FWHM2sigma(add_FWHM);
meanfn = @(x) applyconvfield(x, Sig, add_FWHM, true(Dim), truncation, meanonlat.xvals);
imagesc(meanonlat)

data_info.xvals_vecs = meanonlat.xvals;
data_info.lat_signal = meanonlat;

data_info.npeaks = length(centre_locs);
data_info.maxlocs = cell(1,data_info.npeaks);
for I = 1:data_info.npeaks
    data_info.maxlocs{I} = findlms( meanfn, centre_locs{I} - 1, centre_locs{I}, centre_locs{I} + 1 );
end
max_val = meanfn(data_info.maxlocs{1});
data_info.maxlocs = cell2mat(data_info.maxlocs);

% Scale to ensure that the snr is 1
data_info.lat_signal = meanonlat*(1/max_val);
data_info.meanfn = @(x) meanfn(x)/max_val;
data_info.meanonfinelat = meanonfinelat*(1/max_val);

data_info.lat_std = 1;
data_info.dim = Dim;
data_info.D = length(centre_locs{1});
data_info.FWHM = FWHM;
data_info.nsubj = nsubj;

end

