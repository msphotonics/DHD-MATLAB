function est = align_global_phase(truth, est, mask)
% ALIGN_GLOBAL_PHASE  Remove the arbitrary global phase for display.
%   EST = ALIGN_GLOBAL_PHASE(TRUTH, EST, MASK) multiplies EST by a
%   unit-magnitude phasor chosen so that the cross-correlation
%   sum(conj(TRUTH) .* EST) over the pupil is real and positive.
%
%   This operation is for DISPLAY PURPOSES ONLY.  It does not affect
%   the field-estimated Strehl ratio S_F (which is phase-invariant by
%   construction) and should not be applied before calling
%   FIELD_STREHL.
%
%   The global phase of a digitally holographic reconstructed field
%   estimate is arbitrary because the reference phase is not
%   calibrated in absolute terms.  Aligning to the truth field here
%   is equivalent to what a downstream wavefront sensor or matched
%   filter would do in practice.
%
%   Inputs:
%       truth - Truth pupil-plane complex field, N x N
%       est   - Reconstructed field estimate (before phase alignment), N x N
%       mask  - Binary pupil mask, N x N
%
%   Output:
%       est   - Phase-aligned reconstructed field estimate, N x N

in  = mask > 0;
c   = sum(conj(truth(in)) .* est(in));
est = est * exp(-1j * angle(c));
end
