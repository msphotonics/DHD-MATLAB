function U_P = fresnel_prop_object_to_pupil(U_O, D, delta, wvl, z_P)
% FRESNEL_PROP_OBJECT_TO_PUPIL  Propagate an object field to the pupil plane.
%
%   U_P = FRESNEL_PROP_OBJECT_TO_PUPIL(U_O, D, DELTA, WVL, Z_P) Fresnel
%   propagates the N-by-N object-plane field U_O a distance Z_P, then
%   applies a thin-lens transmittance exp(-jk/(2*z_P)*(x1^2+y1^2)) to
%   model the collecting lens that collimates the light at the pupil plane.
%
%   This models the two-step process described in the chapter background:
%     1. Free-space Fresnel propagation from object to pupil (distance z_P)
%     2. Thin-lens transmittance at the pupil plane (collimating lens)
%
%   The pupil mask should be applied by the calling script after this
%   function returns.
%
%   Inputs:
%       U_O   - Object-plane complex field, N x N
%       D     - Exit-pupil diameter [m]
%       delta - Grid spacing [m]  (= D/N)
%       wvl   - Wavelength [m]
%       z_P   - Propagation distance [m]
%
%   Output:
%       U_P   - Pupil-plane complex field, N x N  [same units as U_O]
%
%   See also FRESNEL_PROP_PUPIL_TO_IMAGE, FFT2C, IFFT2C.

N     = size(U_O, 1);
kwave = 2*pi/wvl;

% Step 1: Fresnel propagate U_O from object plane to pupil plane
% Zero-pad to M = 2N to avoid aliasing in the Fresnel kernel
M  = 2*N + rem(N, 2);
nn = round(N/2) + (1:N);

U_Og = zeros(M);
U_Og(nn, nn) = U_O;

% Forward FFT of padded field
U_Ogt = delta^2 * fftshift(fft2(fftshift(U_Og)));

% Fresnel propagation kernel (1-D, applied separably)
const  = kwave / (2*z_P) * delta^2;
nvec   = ceil(-M/2) : floor((M-1)/2);
kernel = exp(1j * const * nvec(:).^2);

% FFT of kernel
kernelt = delta * fftshift(fft(fftshift(kernel)));

% Multiply by separable 2-D kernel
U_Ogt = (kernelt * kernelt.') .* U_Ogt;

% Inverse FFT and scale to physical units
U_Pm = 1/(1j*wvl*z_P) * exp(1j*kwave*z_P) ...
       * (1/delta)^2 * ifftshift(ifft2(ifftshift(U_Ogt)));

% Extract N x N result
U_Pm = U_Pm(nn, nn);

% Step 2: Apply thin-lens transmittance at pupil plane (collimating lens)
x1vec = (-N/2 : N/2-1)*delta;
[x1, y1] = meshgrid(x1vec, x1vec);
U_P = exp(-1j*kwave/(2*z_P)*(x1.^2 + y1.^2)) .* U_Pm;
end
