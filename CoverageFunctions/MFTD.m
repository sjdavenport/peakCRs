function [distbn, distbn_nobounds, det_hessian, det_hessian_nobounds] = MFTD( peakderiv2, covmate, N, niters )
% MFTD( peakderiv2, covmate, N, niters ) approximates the
% distribution of the first term in the taylor distribution of the maximum.
% ONLY STATIONARY ATM!!
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  peakderiv2   a D by D matrix of the second derivative
%  covmate
%  N
%  niters 
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D
% Lambda = 0.0779; N = 100; mu = -0.1219; Omega = 0.0049*40;
% Delta = Lambda/4; covmate = [Lambda, Delta; Delta, Omega];
% distbn= MFTD( mu, covmate, N );
% histogram(distbn); hold on; histogram(distbn_indep)
% 
% covmate = [ 0.0779, 0; 0, 0.0049*40]
% N = 100; mu = -0.1219; samples = mu + normrnd(0,Omega,1,N);
% mftd_dist = MFTD( mean(samples), covmate, N );
% orig_dist = normrnd(0,sqrt(covmate(1,1)),1, length(mftd_dist) )/sqrt(N);
% clf
% orig_hist = histogram(orig_dist/mean(samples))
% hold on
% histogram(mftd_dist, 'BinWidth', orig_hist.BinWidth)
% alpha = 0.05;
% (1/2)*(prctile(a/mean(samples), 100*(1-alpha) ) - prctile(a/mean(samples), 100*(alpha) ))
% (1/2)*(prctile(mftd_dist, 100*(1-alpha) ) - prctile(mftd_dist, 100*(alpha) ))
%
% % 
% Lambda = eye(2)*0.0779; mu = -0.1219; Omega = eye(3)*0.0049*40;
% Delta = zeros(2,3); covmate = [Lambda, Delta; Delta', Omega];
% MFTD( eye(2), covmate, N, niters )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
   % default option of opt1
%    niters = 100000000;
    niters = 100000;
%     niters = 1000000*10
end

% Obtain the number of dimensions
D = size(peakderiv2,1);

% Ensure that the covariance matrix has the right size
if size(covmate,1) ~= (D + D*(D+1)/2)
    error('The size of the second derivative vector and the covariance matrix are incompatible')
end

% Vectorize the 2nd derivative and convert to a row vector
peakderiv2_vec = vech(peakderiv2)';

%%  Main Function Loop
%--------------------------------------------------------------------------
% joint = mvnrnd([zeros(1,D), peakderiv2], covmate/N, niters)
% Done to avoid numerical issues
joint = mvnrnd([zeros(1,D), peakderiv2_vec*sqrt(N)], covmate, niters)/sqrt(N);
% To fix the numerical issue could generate the data as mu + Sigma^(1/2)X where 
% X is standard multivariate Gaussian and Sigma^(1/2) = sqrtm(covmate).
det_hessian = NaN;
if D == 1
    num = joint(:,1);
    denom = joint(:,2);
    distbn_nobounds = num./denom;
    num = num(denom < peakderiv2_vec/20);
    denom = denom(denom < peakderiv2_vec/20);
    distbn = num./denom;
else
    % Initialize the ellvalue distribution vector
    distbn_nobounds = zeros(niters, D); 
    det_hessian = zeros(niters,1);
    for I = 1:niters
       % Calculate the hessian
       hessian = vech(joint(I,(D+1):end));
       det_hessian(I) = det(hessian);
       
       % Calculate a value of from the MC distribution of thetahat-theta
       distbn_nobounds(I, :) = (inv(hessian)*joint(I,1:D)')'; 
    end
    distbn = distbn_nobounds(det_hessian > det(peakderiv2)/20, :);
    det_hessian_nobounds = det_hessian;
    det_hessian = det_hessian(det_hessian > det(peakderiv2)/20);
    distbn_nobounds = NaN;
end

end

% b = normrnd(0,1,1,niters);
% 
% denom = (peakderiv2 + b*(covmate(2,2)/N)^(1/2));
% denom = denom(denom < peakderiv2/10);
% 
% a = normrnd(0,1,1,length(denom));
% 
% distbn_indep = (covmate(1,1)^(1/2)/sqrt(N))*a./denom;

