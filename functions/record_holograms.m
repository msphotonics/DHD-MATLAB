function [mH, dH_ideal, dH_noisy] = record_holograms(U_S_fine, U_R_vec, K, ...
    sigma_r, usePoisson, ell, enforceWellDepth, g_AD, nBits)
% RECORD_HOLOGRAMS  Record one or more phase-shifted holograms.
%   [MH, DH_IDEAL, DH_NOISY] = RECORD_HOLOGRAMS(U_S_FINE, U_R_VEC,
%   K, SIGMA_R, USEPOISSON, ELL, ENFORCEWELLDEPTH, G_AD, NBITS)
%   records holograms for each reference field in U_R_VEC.
%
%   Inputs:
%       U_S_FINE  - Signal field on the fine (K*N x K*N) FPA grid [sqrt(pe)]
%       U_R_VEC   - Row vector of reference complex amplitudes [sqrt(pe)].
%                   Length 1 for IPRG/PPRG (single hologram).
%                   Length 4 for PSRG (four phase-shifted holograms).
%                   For the IPRG/PPRG, U_R_VEC is the full fine-grid
%                   reference field (K*N x K*N); pass it as a 1-element
%                   cell: {U_R_fine}.
%                   For the PSRG, U_R_VEC is a 1x4 complex scalar vector
%                   (spatially uniform reference).
%       K         - Subpixels per FPA pixel (fine-to-coarse ratio)
%       sigma_r   - Read-noise std dev [pe]
%       usePoisson- true -> Poisson shot noise; false -> Gaussian approx
%       ell       - Pixel well depth [pe]
%       enforceWellDepth - true -> clip at ell before read noise
%       g_AD      - A/D gain [DN/pe]
%       nBits     - A/D bit depth
%
%   Outputs:
%       mH       - Mean hologram irradiance [pe], N x N x nMeas
%       dH_ideal - Noiseless (quantized only) hologram [DN], N x N x nMeas
%       dH_noisy - Noisy hologram [DN], N x N x nMeas
%
%   For single-measurement geometries (IPRG, PPRG), pass U_R_VEC as
%   a cell array containing the full spatial reference field:
%       [mH, dH_ideal, dH_noisy] = record_holograms(U_S_fine, {U_R_fine}, ...)
%   The outputs are then N x N x 1; index as mH(:,:,1) or squeeze(mH).
%
%   For the PSRG, pass a 1x4 scalar complex vector:
%       A_R = sqrt(m_R);
%       U_R_vec = A_R * exp(-1j * [0, pi/2, pi, 3*pi/2]);
%       [mH, dH_ideal, dH_noisy] = record_holograms(U_S_fine, U_R_vec, ...)
%
%   Reference: Section 6.2 (PSRG) and surrounding discussion.

% Determine the number of measurements and reference type
if iscell(U_R_vec)
    % Full spatial reference field(s): cell array of fine-grid arrays
    nMeas     = numel(U_R_vec);
    spatialRef = true;
else
    % Scalar complex amplitudes (PSRG uniform reference)
    nMeas     = numel(U_R_vec);
    spatialRef = false;
end

% Coarse-pixel grid size
N = size(pixel_bin(abs(U_S_fine).^2, K), 1);

mH       = zeros(N, N, nMeas);
dH_ideal = zeros(N, N, nMeas);
dH_noisy = zeros(N, N, nMeas);

for k = 1:nMeas
    if spatialRef
        U_R_k = U_R_vec{k};
    else
        % Scalar amplitude applied uniformly across the fine grid
        U_R_k = U_R_vec(k) * ones(size(U_S_fine));
    end

    mH_k         = pixel_bin(abs(U_S_fine + U_R_k).^2, K);
    mH(:,:,k)    = mH_k;
    dH_ideal(:,:,k) = quantize_dn(mH_k, g_AD, nBits);
    dH_noisy(:,:,k) = detect_dn(mH_k, sigma_r, usePoisson, ell, ...
                                 enforceWellDepth, g_AD, nBits);
end
end
