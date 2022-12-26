function [ellval, inell] = inellipse( x, Sigma, mu, k )
% inellipse( x, Sigma, mu, k ) tests if the point x is in the ellipse given
% by {y: (y - mu)^T * Sigma * (y - mu)^T <= k}.
%--------------------------------------------------------------------------
% ARGUMENTS
%   x      - a column vector representing the point to test.
%   Sigma  - a positive definite matrix representing the covariance matrix 
%             of the ellipse.
%   mu     - a column vector representing the mean of the ellipse. If not 
%               provided, the default value is 0.
%   k      - a scalar representing the threshold value. If not provided, 
%               inell is set to NaN.
%--------------------------------------------------------------------------
% OUTPUT
%   ellval - a scalar representing the value of (y - mu)^T * Sigma * (y - mu)^T 
%           for the point x.
%   inell  - a boolean indicating whether the point x is inside the ellipse. 
%           If k is not provided, inell is set to NaN.
%--------------------------------------------------------------------------
% EXAMPLES
% y = [1,2]'
% Sigma = eye(2);
% mu = [0,0]'
% inellipse( y, Sigma, mu )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 3
    mu = 0;
end

ellval = (x - mu)'*Sigma*(x - mu);
if nargin > 3
    if (x - mu)'*Sigma*(x - mu) <= k
        inell = 1;
    else
        inell = 0;
    end
else
    inell = NaN;
end

end

