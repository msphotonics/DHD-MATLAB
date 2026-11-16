function [UIe, F] = reconstruct_estimate_pprg(holo, M, rowIndex, colIndex, mask)
% RECONSTRUCT_ESTIMATE_PPRG  Obtain a reconstructed field estimate from
%   an off-axis PPRG hologram.
%
%   Inputs:
%       holo      - Zero-mean hologram in photoelectrons, N×N
%       M         - Zero-padded Fourier grid size (M = q_P * N)
%       rowIndex  - Row indices of the direct sideband in F (length N)
%       colIndex  - Col indices of the direct sideband in F (length N)
%       mask      - Binary pupil mask, N×N (unused here, kept for API)
%
%   Outputs:
%       UIe  - Windowed sideband (image-plane estimate), N×N
%       F    - Full M×M Fourier-plane array (for display)

% M×M zero-padded FFT
N   = size(holo, 1);
hp  = zeros(M, 'like', holo);
ie  = centered_indices(M, N);
hp(ie, ie) = holo;
F   = fft2c(hp);

% Extract N×N sideband using rowIndex/colIndex
UIe = F(rowIndex, colIndex);
end
