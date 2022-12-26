function [sig_store, meanfn] = peakgen1D(x, peakspec, peakparams, peakheights, smo)
% PEAKGEN( x, peakspec, aparams, bparams, smo ) generates 1D data with
% peaks from the beta distribution
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  x    the vector at which to evaluate the signal
%  peakspec   a cell array each entry being a 2D vector with start and
%             ending vertices for the peak to occur
%  peakparams  a cell array each entry being a 2D vector specifying the
%              peak parameters (used to generate the peak from the beta pdf)
%  peakheights  a vector giving the height of each peak
% Optional
%  smo   the amount of additional smoothing to apply to the data  (not
%        implemented yet!)
%--------------------------------------------------------------------------
% OUTPUT
%   sig_store   a vector with the same length as x that contains the values
%               of the signal
%--------------------------------------------------------------------------
% EXAMPLES
% peakspec = {[1,4], [5,10], [15,16]}; peakparams = {[1.5,2], [2,2], [2,1.5]}
% smo = 3;
% x = 0:0.1:20; sigstore = peakgen1D( x, peakspec, peakparams, 1, smo );
% plot(x, sigstore)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
%% Set Optional inputs
if ~exist('peakheights', 'var')
    peakheights = 1;
end

if ~exist('smo', 'var')
    smo = 1;
end

% Obtain the increment
increm = x(2) - x(1);

smo = smo/increm;

%%  Error checking
%--------------------------------------------------------------------------
npeaks = length(peakspec);
if length(peakparams) == 1
    peakparams = repmat(peakparams, 1, npeaks);
end

if npeaks ~= length(peakparams)
    error('The number of entries in peakspec and peakparams should be the same')
end

if length(peakheights) == 1
    peakheights = repmat(peakheights, 1, npeaks);
elseif length(peakheights) ~= npeaks
    error('The number of entries in peakspec and peakheights should be the same')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
sig_store = zeros(1, length(x));

meanfn = @(x) 0;
for I = 1:length(peakspec)
    % Fit the linear regression that sends the start of the peak to 0 and
    % the end of the peak to 1. 
    % I.e. solves alpha*x_1 + beta = 0 and alpha*x_2 + beta = 0 
    coeffs201 = inv([1,peakspec{I}(1);1,peakspec{I}(2)])*[0,1]';
    sig_store = sig_store + peakheights(I)*betapdf(coeffs201(2)*x + coeffs201(1), peakparams{I}(1), peakparams{I}(2)); 
    meanfn = @(x) meanfn(x) + peakheights(I)*betapdf(coeffs201(2)*x + coeffs201(1), peakparams{I}(1), peakparams{I}(2));
end

sig_store = fconv(sig_store, smo);

end

