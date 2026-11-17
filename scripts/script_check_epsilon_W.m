% script_check_epsilon_W.m
%
% Encircled-energy check for Section 7 of "Digital-Holographic Detection."
%
% For the bulls-eye (extended) object, Fig. 20(d) shows a constant offset
% in Delta S/N of about -3.1% (FPA in the image plane) and -4.8% (FPA in
% the pupil plane). This script tests whether the offset is the fraction
% of the signal energy that falls outside the image-plane support of
% diameter W, on which the numerical S/N is evaluated:
%
%   epsilon_W = (image energy within W) / (signal energy collected by the FPA)
%
%   FPA in the image plane (IPRG, IP-PSRG): the FPA collects only the image
%     energy that falls on its N x N pixels.
%   FPA in the pupil plane (PPRG, PP-PSRG): the FPA collects all of the
%     energy that passes through the exit pupil, including the image
%     energy outside the W x W region recovered by the filter or by
%     zero padding.
%
% The parameters, aberrations, and bulls-eye object match script_comparison.m.
% The FPA integrates the field over each pixel, so the pixel-integrated
% values are the ones to compare. Expected output (100 realizations):
% 1 - epsilon_W of about 3.2% (FPA in the image plane) and 4.8% (FPA in the
% pupil plane); the field-sample values are printed for reference.

clear; close all; clc;
addpath('../functions');
if exist('rng','file'), rng(1); else, rand('seed',1); randn('seed',1); end %#ok<RAND>

% ------------------------------------------------------------
% Parameters (match script_comparison.m)
% ------------------------------------------------------------
N = 256;  K = 2;  q_I = 3;  q_P = 3;
D = 0.30;  wvl = 1.0e-6;  phi_rms = 1.0;  phaseCorr = D/10;
nReals = 100;

Nfine   = K*N;
delta   = D/N;  p = delta;  W = D;  delta_f = delta/K;
z_I     = q_I*D*p/wvl;
z_P     = q_P*W*p/wvl; %#ok<NASGU>  (equal to z_I here; unit magnification)

xf = (-Nfine/2:Nfine/2-1)*delta_f;
[Xf,Yf] = meshgrid(xf,xf);
pupilMask_fine = double(hypot(Xf,Yf) <= D/2);
pupilMask      = pupilMask_fine(1:K:end,1:K:end);
supp           = pupilMask > 0;          % image support W (diameter N px)
r = hypot(Xf,Yf);
O_fine = double(r <= D/2) - double(r <= 3*D/8) + double(r <= D/4) - double(r <= D/8);

Npad = q_P*N;
ip   = centered_indices(Npad, N);

% ------------------------------------------------------------
% Monte Carlo over aberrations and speckle
% ------------------------------------------------------------
eps_I_samp = zeros(nReals,1);  eps_I_bin = zeros(nReals,1);
eps_P_samp = zeros(nReals,1);  eps_P_bin = zeros(nReals,1);
frac_crop  = zeros(nReals,1);

for ri = 1:nReals
    phi_fine = make_isoplanatic_phase(Nfine,delta_f,phaseCorr,phi_rms,pupilMask_fine);
    P_fine   = pupilMask_fine.*exp(1j*phi_fine);
    R_O_fine = O_fine.*exp(1j*(-pi+2*pi*rand(Nfine)));

    U_P = P_fine.*fresnel_prop_object_to_pupil(R_O_fine,D,delta_f,wvl,z_I);
    U_I = fresnel_prop_pupil_to_image(U_P,D,delta_f,wvl,z_I);

    % FPA in the image plane: energy within W over energy on the FPA
    T  = U_I(1:K:end,1:K:end);
    Tb = pixel_bin(U_I,K);
    eps_I_samp(ri) = sum(abs(T(supp)).^2)/sum(abs(T(:)).^2);
    eps_I_bin(ri)  = sum(abs(Tb(supp)).^2)/sum(abs(Tb(:)).^2);

    % FPA in the pupil plane: image from the pupil field (zero pad by q_P)
    for b = 1:2
        if b == 1, Tp = U_P(1:K:end,1:K:end); else, Tp = pixel_bin(U_P,K); end
        Up = zeros(Npad);  Up(ip,ip) = Tp;
        I  = fft2c(Up);  Ic = I(ip,ip);
        e  = sum(abs(Ic(supp)).^2)/sum(abs(I(:)).^2);
        if b == 1
            eps_P_samp(ri) = e;
            frac_crop(ri)  = sum(abs(Ic(:)).^2)/sum(abs(I(:)).^2);
        else
            eps_P_bin(ri) = e;
        end
    end
end

% ------------------------------------------------------------
% Report
% ------------------------------------------------------------
pr = @(s,v) fprintf('  %-48s %7.4f   1 - eps = %5.2f%%   (std %.4f)\n', ...
    s, mean(v), 100*(1-mean(v)), std(v));
fprintf('\n epsilon_W check, %d realizations\n', nReals);
pr('FPA in image plane, field samples',             eps_I_samp);
pr('FPA in image plane, pixel-integrated field',    eps_I_bin);
pr('FPA in pupil plane, field samples',             eps_P_samp);
pr('FPA in pupil plane, pixel-integrated field',    eps_P_bin);
pr('FPA in pupil plane, energy in the W x W crop',  frac_crop);
fprintf('\n Compare with the Delta S/N offsets in Fig. 20(d):\n');
fprintf('  about -3.1%% (IPRG, IP-PSRG) and -4.8%% (PPRG, PP-PSRG).\n\n');
