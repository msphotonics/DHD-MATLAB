function SF = field_strehl(truth, est, mask)
% FIELD_STREHL  Field-estimated Strehl ratio.
%   SF = FIELD_STREHL(TRUTH, EST, MASK) computes the field-estimated
%   Strehl ratio S_F between the truth pupil-plane complex field TRUTH
%   and the reconstructed field estimate EST, evaluated over the
%   pupil support defined by MASK.
%
%   The field-estimated Strehl ratio is defined as:
%       S_F = |<U_P . conj(U_hat_P)>|^2 /
%             (<|U_P|^2> . <|U_hat_P|^2>)
%
%   where the angle brackets denote spatial averages over the pupil.
%   This metric is phase-invariant (a global phase offset between
%   TRUTH and EST does not affect S_F) and equals 1 when the
%   reconstructed field estimate is a perfect copy of the truth.
%
%   Note: S_F is a property of the reconstructed field estimate, not
%   of the imaging system alone.  Its expected value is related to
%   the signal-to-noise ratio by Eq. (91) of the chapter.
%
%   Inputs:
%       truth  - Truth pupil-plane complex field, N x N
%       est    - Reconstructed field estimate, N x N
%       mask   - Binary pupil mask (1 inside pupil, 0 outside), N x N
%
%   Output:
%       SF     - Field-estimated Strehl ratio in [0, 1]
%
%   Reference: Eq. (90)-(91) of the chapter.

in = mask > 0;
U  = truth(in);
V  = est(in);
SF = abs(mean(U .* conj(V)))^2 / (mean(abs(U).^2) * mean(abs(V).^2));
end
