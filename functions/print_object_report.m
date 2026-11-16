function print_object_report(name, r)
% PRINT_OBJECT_REPORT  Standard console summary for one object.
%   PRINT_OBJECT_REPORT(NAME, R) prints the same block for every geometry.
%
%   Fields of R:
%     mS, mS_meas       requested / measured mean signal pe (full-FPA mean)
%     mR, mR_meas       requested / measured mean reference pe
%     holoMaxPe         max noiseless hologram value [pe]
%     holoMaxDN, maxDN  max recorded hologram value and full scale [DN]
%     eta               pixel-aperture efficiency (PIXEL_ETA)
%     noisePred, noiseMeas, nEff   noise variance: closed-form, numerical,
%                       number of samples (sets the expected scatter)
%     sigPred, sigMeas  signal power: closed-form, numerical
%     snrClosed, snrNum closed-form and numerical SNR
%     kappa             SNR scale factor from where SNR is defined to the
%                       pupil support where S_F is evaluated
%     SF0, SF           noiseless and noisy S_F (pupil plane)
%
%   Predicted S_F = 1/(1 + 1/(kappa*SNR)).
%
%   See also PRINT_RUN_HEADER.

fprintf('\n--- %s %s\n', name, repmat('-',1,max(0,64-numel(name))));
fprintf('  m_S    : requested = %.3f pe, measured = %.3f pe (full-FPA mean)\n', r.mS, r.mS_meas);
fprintf('  m_R    : requested = %.1f pe, measured = %.1f pe\n', r.mR, r.mR_meas);
fprintf('  Holo   : max = %.1f pe, %d of %d DN\n', r.holoMaxPe, round(r.holoMaxDN), r.maxDN);
fprintf('  eta    : %.4f\n', r.eta);
fprintf('  Noise  : predicted = %.4e, measured = %.4e (ratio = %.4f, +/-%.1f%% expected)\n', ...
    r.noisePred, r.noiseMeas, r.noiseMeas/r.noisePred, 100/sqrt(r.nEff));
fprintf('  Signal : predicted = %.4e, measured = %.4e (ratio = %.4f)\n', ...
    r.sigPred, r.sigMeas, r.sigMeas/r.sigPred);
fprintf('  SNR    : closed-form = %.3f, numerical = %.3f, error = %+.2f%%\n', ...
    r.snrClosed, r.snrNum, 100*(r.snrNum-r.snrClosed)/r.snrClosed);
fprintf('  S_F    : noiseless = %.4f, noisy = %.4f\n', r.SF0, r.SF);
fprintf('           predicted = %.4f (closed-form SNR), %.4f (numerical SNR), kappa = %.3f\n', ...
    1/(1+1/(r.kappa*r.snrClosed)), 1/(1+1/(r.kappa*r.snrNum)), r.kappa);
end
