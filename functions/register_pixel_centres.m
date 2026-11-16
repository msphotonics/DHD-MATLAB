function d = register_pixel_centres(d, K)
% REGISTER_PIXEL_CENTRES  Put binned detector data on the nominal grid.
%   D = REGISTER_PIXEL_CENTRES(D, K) shifts each N x N frame of D by
%   a = (K-1)/(2K) detector pixels in x and y with a sub-pixel Fourier shift.
%
%   PIXEL_BIN averages the K x K sub-pixels of each detector pixel, so each
%   binned value represents the pixel CENTRE, which lies (K-1)/2 fine samples
%   from the first sub-pixel. The truth is decimated at the first sub-pixel
%   (P_fine(1:K:end,1:K:end)), so without this step the estimate and truth
%   are sampled at different points and the noiseless S_F is biased low.
%
%   The shift is all-pass (unit magnitude in frequency), so it leaves the
%   noise variance and signal power unchanged. The Nyquist bin is left
%   unshifted so real data stay real.
%
%   See also PIXEL_BIN.

if K == 1, return; end
[N1,N2,nF] = size(d);
a = (K-1)/(2*K);
H = shift_phase(N1,a)*shift_phase(N2,a).';
for f = 1:nF
    d(:,:,f) = real(ifft2(fft2(d(:,:,f)).*H));
end
end

function p = shift_phase(n,a)
% Phase ramp for x[m] -> x[m-a] in standard FFT ordering.
k = 0:n-1;
k(k >= n/2) = k(k >= n/2) - n;
p = exp(-1j*2*pi*k*a/n);
if mod(n,2) == 0, p(n/2+1) = 1; end
p = p(:);
end
