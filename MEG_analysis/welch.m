function [ power_spectrum, ps_fields, real_fft_fields, imag_fft_fields ] = ...
    welch( time_series, segment_length, sample_freq, resadd, scale )
% WELCH() implements Welch's method on a time series
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% time_series
% segment_length
% sample_freq
% resadd
% scale 
% Optional
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
% Set the default frequency
if ~exist('sample_freq', 'var')
    sample_freq = 1;
end

if ~exist('resadd', 'var')
    resadd = 0;
end

if ~exist('scale', 'var')
    scale = 1;
end

% Get the distance between time points.
dt = 1/sample_freq;

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the number of frequencies based on the input segment length
% The fft is of length segment_length but it repeats so you only need to
% take the first half of the frequencies
nfreqs = floor(segment_length/2);

% Obtain the total Time
T = segment_length/sample_freq;

v = 0.05; % The value of the window at the edges

% Calculate the parameter of the Gaussian window so that it gives a
% reasonable covering of the space.
G_window_sigma2 = (-((T-dt)/2)^2)/log(v);

% Calculate the sigma and FWHM for the smoothing kernel in frequency space
sigma_smoothing = sqrt(1/2/G_window_sigma2)/pi;
FWHM = sigma2FWHM(sigma_smoothing);

% Set the parameters for performing the smoothing
params = ConvFieldParams( FWHM, resadd, 0 );

freqs = ((0:(nfreqs-1))/segment_length)*sample_freq;
freq_gap = sample_freq/segment_length; % The gap between neighbouring frequencies
max_freq = max(freqs);

% % Initialize the fields for the fourier transform
% real_fft_fields = Field(true(nfreqs,1));
% imag_fft_fields = Field(true(nfreqs,1));

% Calculate the number of fields
nfields = floor(length(time_series)/nfreqs) - 1;

% Take into account the edge effect
extra = ceil(8*FWHM2sigma(FWHM));

% Initialize the fields for the real and imaging ffts
real_fft_fields = Field(zeros(nfreqs + 2*extra, nfields), true(nfreqs+ 2*extra,1));
imag_fft_fields = Field(zeros(nfreqs + 2*extra, nfields), true(nfreqs+ 2*extra,1));
real_fft_fields.xvals{1} = (-extra*freq_gap):freq_gap:(max_freq+extra*freq_gap);
% imag_fft_fields.xvals{1} = freqs;

%%  Main Function Loop
%--------------------------------------------------------------------------
start = 1;
%
% if window == 1
%     window = hamming(segment_length)';
% else
%     window = 1;
% end

% Ensure that no enlarging takes place (for now at least)
params.enlarge = 0;

% Loop to calculate the ffts of the individual segments
for I = 1:nfields
    % Calculate the fft of the Ith segment
    %     segment_fft = fft(time_series(start:start+segment_length-1).*window);
    segment_fft = fft(time_series(start:start+segment_length-1));
    
    real_fft = real(segment_fft(1:nfreqs));
    imag_fft = imag(segment_fft(1:nfreqs));
    
    % Separate into real and imaginary parts
    real_fft_fields.field(:,I) = [fliplr(real_fft(1:extra)), real_fft, fliplr(real_fft((end-extra+1):end))];
    imag_fft_fields.field(:,I) = [fliplr(imag_fft(1:extra)), imag_fft, fliplr(imag_fft((end-extra+1):end))];
    
    % Update the start point
    start = start + floor(segment_length/2);
end

% Smooth the fields and cutoff the bit that was added to correct for edge
% effects
real_smoothed_fields = cut_field(convfield(real_fft_fields, params), extra)*(1/T);
imag_smoothed_fields = cut_field(convfield(imag_fft_fields, params), extra)*(1/T);

% Scale sum the sum of the kernel to ensure everything is even
if scale == 1
    fieldofones = Field(ones(nfreqs + 2*extra, 1), true(nfreqs + 2*extra, 1));
    fieldofones.xvals{1} = (-extra*freq_gap):freq_gap:(max_freq+extra*freq_gap);
    convolvedones = cut_field(convfield(fieldofones, params), extra);
    real_smoothed_fields.field = real_smoothed_fields.field./convolvedones.field;
    imag_smoothed_fields.field = imag_smoothed_fields.field./convolvedones.field;
end

% Calculate the individual power spectrum fields
ps_fields = real_smoothed_fields.^2 + imag_smoothed_fields.^2;
ps_fields.xvals = {0:(freq_gap/(resadd+1)):max_freq}; % Need to change cutfield so that this line is unnecessary

% Obtain the overall power spectrum as the mean of the individual ones
power_spectrum = mean(ps_fields);

end

% % Other stuff for detects mean I think, probably deprecated
% real_mean = mean(real_fft_field_data,2);
% imag_mean = mean(imag_fft_field_data,2);
%
% real_demeaned = real_fft_field_data - real_mean;
% imag_demeaned = imag_fft_field_data - imag_mean;
%
% dm_ps_fields = Field(true(nfreqs,1));
% dm_ps_fields.field = sqrt(real_demeaned.^2 + imag_demeaned.^2);
% dm_power_spectrum = mean(dm_ps_fields.field,2);

