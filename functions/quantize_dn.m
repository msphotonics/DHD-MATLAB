function dH = quantize_dn(meanPe, g, nBits)
% QUANTIZE_DN  Noiseless analogue-to-digital conversion.
%   DH = QUANTIZE_DN(MEANPE, G, NBITS) converts the mean
%   photoelectron count MEANPE to integer digital numbers (DN)
%   by applying A/D gain G [DN/pe], rounding to the nearest
%   integer, and clipping to the range [0, 2^NBITS - 1].
%
%   This models a noiseless (ideal) hologram; use DETECT_DN to
%   include shot noise, read noise, and well-depth clipping.
%
%   Inputs:
%       meanPe  - Mean photoelectron array [pe]
%       g       - A/D gain [DN/pe]; typically (2^nBits-1)/ell
%       nBits   - A/D bit depth (integer)
%
%   Output:
%       dH      - Recorded hologram [integer DN]
%
%   See also DETECT_DN.

dH = min(max(round(g * real(meanPe)), 0), 2^nBits - 1);
end
