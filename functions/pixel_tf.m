function H = pixel_tf(x, K)
% PIXEL_TF  Pixel-aperture transfer function (Dirichlet kernel).
%   H = PIXEL_TF(X, K) returns the transfer function of the
%   spatial averaging performed by PIXEL_BIN over K sub-samples
%   per pixel, evaluated at normalized frequency X (cycles/pixel).
%
%   The result is the Dirichlet kernel:
%       H(x) = sin(pi*x) / (K * sin(pi*x/K))
%
%   which reduces to sinc(x) = sin(pi*x)/(pi*x) as K -> infinity,
%   and to 1 for K = 1 (no pixel aperture effect).
%
%   X is typically x_1/D evaluated at the sideband centre, e.g.:
%       pixel_tf(u1/q + 1/4, K)   for the IPRG or PPRG
%       pixel_tf(0, K)             for the PSRG (DC reference)
%
%   Reference: Eq. (38) of the chapter.

H  = ones(size(x));
nz = abs(sin(pi * x / K)) > eps;
H(nz) = sin(pi * x(nz)) ./ (K * sin(pi * x(nz) / K));
end
