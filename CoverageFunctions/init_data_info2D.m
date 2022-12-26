function data_info = init_data_info2D( FWHM, nsubj, Rad, Smo, Dim, centre_locs, Mag )
% INIT_DATA_INFO2D initializes the data_info framework for running coverage
% testing
%--------------------------------------------------------------------------
% ARGUMENTS
% Mag       is a vector of length npeaks giving the magnitude of the 
%           resulting signal at each peak. Note that if Mag is just a real
%           number rather than a vector then all peaks are taken to have
%           magnitude Mag. npeaks is then determined by the length of
%           centre_locs.
% Rad       is a vector of length npeaks giving the radius of the spheroid 
%           signal at each peak prior to being smoothed. Note that if Rad 
%           is just a real number rather than a vector then all peaks are 
%           taken to have radius Rad.
% Smo       is a vector of length npeaks that gives the smoothing applied
%           to the signal at each peak. The smoothing applied is Gaussian 
%           with the same FWHM in each x,y and z directions. Note that if Smo 
%           is just a real number rather than a vector then all peaks are
%           smoothed with FWHM Smo.
% Dim       a vector of length D giving the dimensions of the output image. 
%           For example Dim =[20,50] means that the output image is 20 x 50.
% centre_locs   is a cell array of length npeaks such that the nth entry is
%               a length D vector giving the coordinates of the centre location
%               of the nth peak.
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
truncation = ceil(4*FWHM2sigma(FWHM));
for I = 1:length(centre_locs)
    centre_locs{I} = centre_locs{I} + truncation;
end 
for I = 1:length(Dim)
    Dim(I) = Dim(I) + 2*truncation;
end
Sig = peakgen( Mag, Rad, Smo, Dim, centre_locs ); add_FWHM = 2;
data_info.Sig = Sig;
resadd = 0; params = ConvFieldParams([add_FWHM,add_FWHM], resadd);
meanonlat = convfield(Sig, params);
params.resadd = 1; meanonfinelat = convfield(Sig, params);
add_truncation = ceil(8*FWHM2sigma(add_FWHM));
meanfn = @(x) applyconvfield(x, Sig, add_FWHM, true(Dim), add_truncation, meanonlat.xvals);
imagesc(meanonlat)

data_info.xvals_vecs = meanonlat.xvals;
data_info.lat_signal = meanonlat;

data_info.npeaks = length(centre_locs);
data_info.maxlocs = cell(1,data_info.npeaks);
for I = 1:data_info.npeaks
    data_info.maxlocs{I} = findlms( meanfn, centre_locs{I}, centre_locs{I} - 1, centre_locs{I} + 1 );
end

max_val = meanfn(data_info.maxlocs{1});
data_info.max_val = max_val;
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

data_info.resadd = 51;

% params.resadd = 51;
% meanonfinelat = convfield(Sig, params);
% lms = lmindices(meanonfinelat.field, 3);
% data_info.maxlocs = xvaleval(lms, meanonfinelat.xvals)

end

