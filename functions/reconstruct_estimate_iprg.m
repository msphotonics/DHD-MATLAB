function [UPe, F] = reconstruct_estimate_iprg(holo, M, rowIndex, colIndex, mask)
% RECONSTRUCT_ESTIMATE_IPRG  Obtain a reconstructed field estimate from
%   an off-axis IPRG hologram.
%
%   [UPE, F] = RECONSTRUCT_ESTIMATE_IPRG(HOLO, M, ROWINDEX, COLINDEX,
%   MASK) zero-pads the N-by-N hologram HOLO (in photoelectrons) into
%   an M-by-M grid, applies a centered inverse FFT to obtain the
%   Fourier plane F, then extracts the direct spatial-heterodyne
%   sideband at indices (ROWINDEX, COLINDEX) and windows it by MASK
%   to produce the pupil-plane reconstructed field estimate UPE.
%
%   The inverse FFT convention places the U_S * U_R* cross-term
%   (direct sideband) at positive spatial frequencies when the
%   reference carries exp(+j*2*pi*(x+y)/(4*p)) phase (cf. Eq. 28-30
%   of the chapter).
%
%   Inputs:
%       holo      - Zero-mean hologram in photoelectrons, N x N
%                   (obtained by dividing the integer-DN hologram by
%                   the A/D gain g_AD)
%       M         - Zero-padded Fourier grid size (M = q_I * N)
%       rowIndex  - Row indices of the direct sideband in F
%       colIndex  - Column indices of the direct sideband in F
%       mask      - Binary pupil mask, N x N
%
%   Outputs:
%       UPe  - Pupil-plane reconstructed field estimate, N x N
%       F    - Full Fourier-plane array, M x M (for display and
%              noise estimation)
%
%   Reference: Eq. (34) of the chapter.
%
%   See also FIELD_STREHL, ALIGN_GLOBAL_PHASE.

N  = size(holo, 1);
hp = zeros(M, 'like', holo);
ie = centered_indices(M, N);
hp(ie, ie) = holo;
F   = ifft2c(hp);
UPe = mask .* F(rowIndex, colIndex);
end
