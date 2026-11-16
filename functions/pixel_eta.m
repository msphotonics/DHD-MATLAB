function eta = pixel_eta(U_cross, K)
% PIXEL_ETA  Pixel-aperture power efficiency of a binned hologram cross term.
%   ETA = PIXEL_ETA(U_CROSS, K) returns the fraction of signal power that
%   survives K x K pixel binning, where U_CROSS = U_S .* conj(U_R) is the
%   hologram cross term on the fine grid (the reference amplitude cancels).
%
%   Pixel binning box-filters U_CROSS. With H the pixel transfer function
%   (PIXEL_TF, in cycles per detector pixel) and U~ = FT{U_CROSS},
%
%       eta = sum(|H|^2 |U~|^2) / sum(|U~|^2).
%
%   Using the cross term (not U_S alone) covers every geometry the same
%   way: an off-axis reference moves the spectrum to the carrier frequency,
%   and a spherical reference removes the image-plane quadratic phase.
%
%   See also PIXEL_TF, PIXEL_BIN.

Nf  = size(U_cross,1);
W2  = abs(fft2c(U_cross)).^2;           % cross-term power spectrum
x   = K*((0:Nf-1)-Nf/2)/Nf;             % cycles per detector pixel
h   = pixel_tf(x,K);
H2  = abs(h(:)*h(:).').^2;              % separable 2-D |H|^2
eta = sum(H2(:).*W2(:))/sum(W2(:));
end
