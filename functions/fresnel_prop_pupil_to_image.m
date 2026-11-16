function U_I = fresnel_prop_pupil_to_image(U_P, D, delta, wvl, z_I)
% FRESNEL_PROP_PUPIL_TO_IMAGE  Propagate a pupil field to the image plane.
%
%   U_I = FRESNEL_PROP_PUPIL_TO_IMAGE(U_P, D, DELTA, WVL, Z_I) applies a
%   thin-lens transmittance to the N-by-N pupil field U_P and then Fresnel
%   propagates it a distance Z_I to obtain the image-plane complex field U_I.
%
%   The thin-lens transmittance exp(-jk/(2*z_I)*(x^2+y^2)) cancels the
%   quadratic phase term in the Fresnel integral, yielding a result that is
%   analogous to Fraunhofer diffraction but computed with higher accuracy via
%   the full Fresnel machinery (cf. Voelz, "Computational Fourier Optics").
%
%   The output U_I is on the SAME N-by-N grid with the SAME spacing DELTA as
%   the input U_P.  No pixel binning is applied; the sampling quotient q_I
%   (set via Z_I) determines how many pixels span the diffraction-limited PSF.
%
%   Inputs:
%       U_P   - Pupil-plane complex field, N x N
%       D     - Exit-pupil diameter [m]
%       delta - Grid spacing [m]  (= D/N)
%       wvl   - Wavelength [m]
%       z_I   - Propagation distance [m]
%
%   Output:
%       U_I   - Image-plane complex field, N x N  [same units as U_P]
%
%   See also FFT2C, IFFT2C.

N     = size(U_P, 1);
kwave = 2*pi/wvl;

% Pupil-plane coordinate grid
x1vec = (-N/2 : N/2-1)*delta;
[x1, y1] = meshgrid(x1vec, x1vec);

% Apply thin-lens transmittance (positive-lens focusing phase)
U_Pp = exp(-1j*kwave/(2*z_I)*(x1.^2 + y1.^2)) .* U_P;

% Zero-pad to M = 2N (guard-band to avoid aliasing in Fresnel kernel)
M  = 2*N + rem(N, 2);
nn = round(N/2) + (1:N);   % embed indices for N points in M-point grid

U_Ppg = zeros(M);
U_Ppg(nn, nn) = U_Pp;

% Forward FFT of padded field
U_Ppgt = delta^2 * fftshift(fft2(fftshift(U_Ppg)));

% Fresnel propagation kernel (1-D, applied separably)
const  = kwave / (2*z_I) * delta^2;
nvec   = ceil(-M/2) : floor((M-1)/2);
kernel = exp(1j * const * nvec(:).^2);

% FFT of kernel
kernelt = delta * fftshift(fft(fftshift(kernel)));

% Multiply by separable 2-D kernel
U_Ppgt = (kernelt * kernelt.') .* U_Ppgt;

% Inverse FFT and scale to physical units
U_Ip = 1/(1j*wvl*z_I) * exp(1j*kwave*z_I) ...
       * (1/delta)^2 * ifftshift(ifft2(ifftshift(U_Ppgt)));

% Extract N x N result (center of M x M grid)
U_I = U_Ip(nn, nn);
end
