function USe = reconstruct_estimate_psrg(dH0_pe, mask)
% RECONSTRUCT_ESTIMATE_PSRG  Obtain a reconstructed signal-field estimate
%   from four phase-shifted PSRG holograms.
%
%   USE = RECONSTRUCT_ESTIMATE_PSRG(DH0_PE, MASK) computes the
%   signal-field estimate at the FPA from the four zero-mean phase-shifted
%   holograms in DH0_PE using the four-step phase-shift algorithm
%   (cf. Eq. 77 of the chapter):
%
%       U_hat_S = C_PS * [ (d1 - d3) + j*(d4 - d2) ]
%
%   where d1..d4 are the holograms for delta = 0, pi/2, pi, 3*pi/2
%   and C_PS = 1/4.
%
%   The result is the signal-field estimate at the FPA — which may be in
%   the pupil plane (PP-PSRG) or the image plane (IP-PSRG). The calling
%   script is responsible for any subsequent propagation to the pupil plane.
%   The mask is applied to suppress out-of-support noise.
%
%   Inputs:
%       dH0_pe  - Zero-mean hologram stack in photoelectrons,
%                 N x N x 4 (slices: delta = 0, pi/2, pi, 3*pi/2)
%       mask    - Binary support mask, N x N (pupil mask for PP-PSRG;
%                 ones(N) for IP-PSRG since support is not known a priori)
%
%   Output:
%       USe  - Signal-field estimate at the FPA, N x N
%
%   Reference: Eq. (77) of the chapter.

h1 = dH0_pe(:,:,1);    % delta = 0
h2 = dH0_pe(:,:,2);    % delta = pi/2
h3 = dH0_pe(:,:,3);    % delta = pi
h4 = dH0_pe(:,:,4);    % delta = 3*pi/2

USe = (1/4) * mask .* ( (h1 - h3) + 1j*(h4 - h2) );
end
