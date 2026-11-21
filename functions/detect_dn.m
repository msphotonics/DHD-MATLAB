function dH = detect_dn(meanPe, sigma_r, usePoisson, ell, enforceWellDepth, g, nBits)
% DETECT_DN  Simulate noisy hologram detection.
%   DH = DETECT_DN(MEANPE, SIGMA_R, USEPOISSON, ELL,
%                  ENFORCEWELLDEPTH, G, NBITS) simulates the
%   complete detection chain:
%       1. Photon noise (Poisson or Gaussian approximation)
%       2. Well-depth clipping (if ENFORCEWELLDEPTH is true)
%       3. Read noise (zero-mean Gaussian, std = SIGMA_R [pe])
%       4. A/D conversion: gain G [DN/pe], rounding (quantization
%          noise variance = 1/12 DN^2), clip to [0, 2^NBITS-1]
%
%   Inputs:
%       meanPe           - Mean photoelectron array [pe]
%       sigma_r          - Read-noise standard deviation [pe]
%       usePoisson       - true  -> Poisson photon noise (requires
%                                   Statistics Toolbox)
%                          false -> Gaussian approximation
%       ell              - Pixel well depth [pe]
%       enforceWellDepth - true  -> clip photoelectrons at ell
%                                   before read noise
%       g                - A/D gain [DN/pe]
%       nBits            - A/D bit depth (integer)
%
%   Output:
%       dH  - Noisy recorded hologram [integer DN]
%
%   Reference: Eq. (22) and surrounding discussion in the chapter.
%
%   See also QUANTIZE_DN.

meanPe = max(real(meanPe), 0);         % ensure non-negative mean

% Step 1: photon noise
if usePoisson
    photon = poissrnd(meanPe);
else
    photon = meanPe + sqrt(meanPe) .* randn(size(meanPe));
end

% Step 2: well-depth clipping
if enforceWellDepth
    photon = min(photon, ell);
end

% Step 3: read noise
analog = photon + sigma_r * randn(size(meanPe));

% Step 4: A/D conversion with rounding and clipping
dH = min(max(round(g * analog), 0), 2^nBits - 1);
end
