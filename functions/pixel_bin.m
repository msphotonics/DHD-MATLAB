function low = pixel_bin(high, K)
% PIXEL_BIN  Integrate a fine-grid irradiance over K-by-K pixel bins.
%   LOW = PIXEL_BIN(HIGH, K) averages each non-overlapping K-by-K
%   block of the fine-grid array HIGH to produce the coarse-pixel
%   array LOW.  This models the spatial integration of photoelectrons
%   within a square detector pixel of pitch p, where the fine grid
%   has pitch p/K.
%
%   HIGH must have dimensions that are integer multiples of K.
%   LOW has dimensions size(HIGH)/K.
%
%   For K = 1 the function returns HIGH unchanged.
%
%   Reference: Section 4.2 of the chapter (pixel-aperture MTF).

[nr, nc] = size(high);
low = squeeze(mean(mean(reshape(high, K, nr/K, K, nc/K), 1), 3));
end
