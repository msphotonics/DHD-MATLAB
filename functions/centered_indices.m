function ind = centered_indices(Nlarge, Nsmall)
% CENTERED_INDICES  Row/column indices of a centered sub-array.
%   IND = CENTERED_INDICES(NLARGE, NSMALL) returns the NSMALL
%   integer indices into a 1-D array of length NLARGE that
%   correspond to the centered sub-array of length NSMALL.
%
%   Useful for zero-padding (embed an N-point array into the
%   center of an M-point array) and cropping (extract the central
%   N-point region from an M-point array).
%
%   Example:
%       A = zeros(768, 768);
%       ip = centered_indices(768, 256);
%       A(ip, ip) = B;   % embed 256x256 B into the center of A

first = floor((Nlarge - Nsmall) / 2) + 1;
ind   = first : first + Nsmall - 1;
end
