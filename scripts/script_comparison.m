% ============================================================
% script_comparison.m
%
% Section 7: Performance comparison of the four recording geometries
%   IPRG     off-axis image-plane recording geometry
%   PPRG     off-axis pupil-plane recording geometry
%   IP-PSRG  on-axis image-plane phase-shifting recording geometry
%   PP-PSRG  on-axis pupil-plane phase-shifting recording geometry
%
% Monte Carlo sweep over the closed-form SNR. For each realization, a new
% aberration (phase screen), a new bulls-eye speckle realization, and new
% detection noise are drawn. For each target SNR, m_S is set so that the
% closed-form SNR of that geometry equals the target exactly.
%
% Produces one 2x2 figure. Each row is measured in one plane and plotted
% against the S/N on the support where it is measured, so the predicted
% curves of all four geometries collapse onto one curve, 1/(1+1/x):
%   Point source, pupil plane, vs pupil-support S/N = kappa_P*S/N
%     (a) field-estimated Strehl ratio S_F
%     (b) Delta S/N (%): measured pupil-support S/N vs kappa_P*S/N
%   Extended (bulls-eye) object, image plane, vs image-support S/N = kappa_I*S/N
%     (c) field-estimated Strehl ratio S_F of the speckled image field
%     (d) Delta S/N (%): measured image-support S/N vs kappa_I*S/N
%
% kappa_P and kappa_I convert each geometry's closed-form S/N to the S/N on
% the support where the metric is measured. Each follows from at most two
% effects:
%   signal concentration  m_S is a full N x N FPA mean, but the signal fills
%                         only the circular support, pi/4 of the FPA -> 4/pi
%   noise filtering       reaching the support keeps only a fraction of the
%                         noise band -> 1/(that fraction)
% The off-axis closed forms already contain their sideband filter (4q^2/pi).
%
% kappa_I assumes all of the image energy lies inside the support W. The
% aberrated PSF blurs roughly 3-5% of the bulls-eye energy past the edge of
% W, so panel (d) shows a constant negative offset for every geometry. This
% loss could be folded into eta as an encircled-energy factor
% eps_W = (energy of the true image inside W)/(total energy), computed from
% the noiseless field like PIXEL_ETA; (d) would then sit near zero. It is
% left out here so that (d) shows the loss explicitly.
%
%            kappa_P (pupil)  kappa_I (image)        reason
%   IPRG     1                4/pi                   image: concentration
%   PPRG     1                1                      Eq. (67) already averages over W
%   IP-PSRG  4q_I^2/pi        4/pi                   pupil: pupil mask keeps pi/(4q_I^2)
%                                                    of the image-plane noise band
%   PP-PSRG  4/pi             (4/pi)(4q_P^2/pi)      pupil: concentration; image: the
%                                                    pupil-to-image transform spreads
%                                                    white pupil noise over q_P*W
%
% The physics, sampling, and processing match the four illustrative
% scripts exactly (fine grid, Fresnel propagation, K x K binning, references,
% REGISTER_PIXEL_CENTRES, PIXEL_ETA, full-FPA m_S, reference-only noise).
%
% All functions are in ../functions/.
% Figures and results are saved to ../outputs/.
% ============================================================

clear; close all; clc;

% Run from this script's folder, so ../functions and ../outputs resolve
% wherever MATLAB's current folder happens to be
cd(fileparts(mfilename('fullpath')));
addpath('../functions');

% ------------------------------------------------------------
% Parameters (match the illustrative scripts)
% ------------------------------------------------------------

N = 256;  K = 2;  q_I = 3;  q_P = 3;  S = 1;

D = 0.30;          wvl = 1.0e-6;
phi_rms = 1.0;     phaseCorr = D/10;
sigma_r = 100;     ell = 100e3;      m_R = 0.75*ell;
nBits = 14;        g_AD = (2^nBits-1)/ell;  sigma_q2 = 1/12;

enforceWellDepth = true;
usePoisson = license('test','statistics_toolbox');

% Monte Carlo
SNR_vec = logspace(-1, 2, 10);   % target closed-form SNR
nReals  = 100;                   % realizations per SNR point

figW = 6.0;  figH = 6.0;  fontSize = 10;   % 6 in x 6 in, as Figs. 2 and 3
exportFigures = true;  saveResults = true;

set(groot,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultColorbarTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',fontSize,'defaultTextFontSize',fontSize, ...
    'defaultColorbarFontSize',fontSize,'defaultAxesLineWidth',0.5);

% ------------------------------------------------------------
% Derived quantities
% ------------------------------------------------------------

Nfine   = K*N;
delta   = D/N;  p = delta;  W = D;  delta_f = delta/K;
z_I     = q_I*D*p/wvl;               % IPRG, IP-PSRG
z_P     = q_P*W*p/wvl;               % PPRG, PP-PSRG (Eq. 53)
kwave   = 2*pi/wvl;
sigma_n2 = m_R + sigma_r^2 + sigma_q2/g_AD^2;

% Fine grid (pupil and image planes share the same sampling)
xf = (-Nfine/2:Nfine/2-1)*delta_f;
[Xf,Yf] = meshgrid(xf,xf);
pupilMask_fine = double(hypot(Xf,Yf) <= D/2);
pupilMask      = pupilMask_fine(1:K:end,1:K:end);
supp           = pupilMask > 0;      % pupil support (diameter D) and image
                                     % support W (diameter W = D): both N px
O_fine         = make_bullseye_target_phys(Xf,Yf,D);

% Image-plane quadratic phase removed by the reconstruction (spherical
% reference, or image formed by transforming the pupil estimate)
Q_I = exp(-1j*kwave/(2*z_I)*(Xf.^2+Yf.^2));
Q_P = exp(-1j*kwave/(2*z_P)*(Xf.^2+Yf.^2));

% ------------------------------------------------------------
% Geometries (references identical to the illustrative scripts)
% ------------------------------------------------------------

names = {'IPRG','PPRG','IP-PSRG','PP-PSRG'};
nG    = numel(names);
G     = cell(1,nG);

% IPRG: spherical wave + tilt, Eq. (28)
x_R = D*q_I/4;
U_R = (1/1j)*exp(1j*kwave*z_I) .* exp(1j*kwave/(2*z_I)*(Xf.^2+Yf.^2)) ...
    .* exp(1j*2*pi*x_R*(Xf+Yf)/(wvl*z_I));
G{1} = offaxis_geometry('IPRG', sqrt(m_R)*U_R/sqrt(mean(abs(U_R(:)).^2)), q_I, N);

% PPRG: tilted plane wave, Eq. (49)
x_R = D*q_P/4;
U_R = exp(-1j*2*pi*x_R*(Xf+Yf)/(wvl*z_P));
G{2} = offaxis_geometry('PPRG', sqrt(m_R)*U_R/sqrt(mean(abs(U_R(:)).^2)), q_P, N);

% IP-PSRG: on-axis spherical wave, four piston shifts
U_R = (1/1j)*exp(1j*kwave*z_I) .* exp(1j*kwave/(2*z_I)*(Xf.^2+Yf.^2));
G{3} = psrg_geometry('IP-PSRG', sqrt(m_R)*U_R/sqrt(mean(abs(U_R(:)).^2)), q_I, N);

% PP-PSRG: on-axis uniform plane wave, four piston shifts
G{4} = psrg_geometry('PP-PSRG', sqrt(m_R)*ones(Nfine), q_P, N);

% Closed-form SNR prefactor: SNR = eta * alpha * m_R * m_S / sigma_n^2
alpha   = [4*q_I^2/pi, 4*q_P^2/pi, 4/S, 4/S];
% Closed-form S/N -> S/N on the support where each metric is measured
kappa_P   = [1, 1, 4*q_I^2/pi, 4/pi];                 % pupil support
kappa_I = [4/pi, 1, 4/pi, (4/pi)*(4*q_P^2/pi)];     % image support W

fprintf('\n%s\n SECTION 7: PERFORMANCE COMPARISON\n%s\n', repmat('=',1,70), repmat('=',1,70));
fprintf('  N = %d, K = %d, q_I = %.3f, q_P = %.3f, S = %d\n', N, K, q_I, q_P, S);
fprintf('  D = %.3f m, lambda = %.2e m, z_I = %.4f m, z_P = %.4f m\n', D, wvl, z_I, z_P);
fprintf('  %d-bit A/D, g_AD = %.4f DN/pe, sigma_n^2 = %.2f pe^2\n', nBits, g_AD, sigma_n2);
fprintf('  %d SNR points from %.2g to %.2g, %d realizations each\n', ...
    numel(SNR_vec), SNR_vec(1), SNR_vec(end), nReals);
fprintf('  kappa_P (pupil) = %.3f (IPRG), %.3f (PPRG), %.3f (IP-PSRG), %.3f (PP-PSRG)\n', kappa_P);
fprintf('  kappa_I (image) = %.3f (IPRG), %.3f (PPRG), %.3f (IP-PSRG), %.3f (PP-PSRG)\n\n', kappa_I);

% ------------------------------------------------------------
% Monte Carlo sweep
% ------------------------------------------------------------

nSNR = numel(SNR_vec);
% Arrays: (geometry, SNR point, realization)
SF_pt   = zeros(nG,nSNR,nReals);   % point source, pupil plane
SNR_pt  = zeros(nG,nSNR,nReals);   % measured pupil-support S/N
SF_be   = zeros(nG,nSNR,nReals);   % bulls-eye, image plane
SNR_be  = zeros(nG,nSNR,nReals);   % measured image-support S/N
sat_pt  = zeros(nG,nSNR,nReals);   % fraction of saturated hologram pixels
sat_be  = zeros(nG,nSNR,nReals);
eta_all = zeros(nG,nReals,2);

tic;
for ri = 1:nReals
    rng(ri);

    % New aberration and speckle realization
    phi_fine = make_isoplanatic_phase(Nfine,delta_f,phaseCorr,phi_rms,pupilMask_fine);
    P_fine   = pupilMask_fine.*exp(1j*phi_fine);
    R_O_fine = O_fine.*exp(1j*(-pi+2*pi*rand(Nfine)));

    % Pupil fields: object -> pupil (collimate), aberrated pupil
    U_Pbe_I = P_fine.*fresnel_prop_object_to_pupil(R_O_fine,D,delta_f,wvl,z_I);
    if z_P == z_I
        U_Pbe_P = U_Pbe_I;
    else
        U_Pbe_P = P_fine.*fresnel_prop_object_to_pupil(R_O_fine,D,delta_f,wvl,z_P);
    end

    % Image fields: pupil -> image (focus)
    U_Ipt   = fresnel_prop_pupil_to_image(P_fine, D,delta_f,wvl,z_I);
    U_Ibe_I = fresnel_prop_pupil_to_image(U_Pbe_I,D,delta_f,wvl,z_I);
    U_Ibe_P = fresnel_prop_pupil_to_image(U_Pbe_P,D,delta_f,wvl,z_P);

    % Truths (alias-free decimation; the estimates are registered to this grid)
    T_pupil = P_fine(1:K:end,1:K:end);              % point source, pupil plane
    T_img_I = U_Ibe_I.*Q_I;  T_img_I = T_img_I(1:K:end,1:K:end);
    T_img_P = U_Ibe_P.*Q_P;  T_img_P = T_img_P(1:K:end,1:K:end);
    T_img   = {T_img_I, T_img_P, T_img_I, T_img_P}; % bulls-eye, image plane

    % Signal field at each FPA: {point source, bulls-eye}
    U_FPA = { U_Ipt,  U_Ibe_I;   % IPRG
              P_fine, U_Pbe_P;   % PPRG
              U_Ipt,  U_Ibe_I;   % IP-PSRG
              P_fine, U_Pbe_P }; % PP-PSRG

    for gi = 1:nG
        g = G{gi};

        % Reference-only recording (signal = 0), carried to both planes
        dN  = record_registered(g, zeros(Nfine), K, sigma_r, usePoisson, ...
            ell, enforceWellDepth, g_AD, nBits);
        nPu = pupil_estimate(g, dN, g_AD, pupilMask);  nPu = nPu(supp);
        nIm = image_estimate(g, dN, g_AD, pupilMask);  nIm = nIm(supp);

        for oi = 1:2
            U   = U_FPA{gi,oi};
            m0  = mean(mean(pixel_bin(abs(U).^2,K)));      % full-FPA mean
            eta = pixel_eta(U.*conj(g.refs{1}),K);
            eta_all(gi,ri,oi) = eta;

            for si = 1:nSNR
                % m_S that gives the target closed-form SNR
                mS  = SNR_vec(si)*sigma_n2/(eta*alpha(gi)*m_R);
                [dH0, sat] = record_registered(g, sqrt(mS/m0)*U, K, sigma_r, ...
                    usePoisson, ell, enforceWellDepth, g_AD, nBits);

                if oi == 1   % point source: pupil plane
                    E = pupil_estimate(g, dH0, g_AD, pupilMask);
                    SF_pt(gi,si,ri)  = field_strehl(T_pupil, E, pupilMask);
                    SNR_pt(gi,si,ri) = support_snr(E(supp), nPu);
                    sat_pt(gi,si,ri) = sat;
                else         % bulls-eye: image plane
                    E = image_estimate(g, dH0, g_AD, pupilMask);
                    SF_be(gi,si,ri)  = field_strehl(T_img{gi}, E, pupilMask);
                    SNR_be(gi,si,ri) = support_snr(E(supp), nIm);
                    sat_be(gi,si,ri) = sat;
                end
            end
        end
    end

    if ri == 1 || mod(ri,10) == 0 || ri == nReals
        fprintf('  realization %3d of %d done (%.0f s elapsed)\n', ri, nReals, toc);
    end
end

% ------------------------------------------------------------
% Statistics
% ------------------------------------------------------------

x_P = kappa_P(:)*SNR_vec;      % predicted pupil-support S/N, nG x nSNR
x_I = kappa_I(:)*SNR_vec;    % predicted image-support S/N, nG x nSNR

mn = @(x) mean(x,3);   sd = @(x) std(x,0,3);
pct = @(meas,pred) 100*(meas - pred)./pred;

R.x_P = x_P;   R.x_I = x_I;
R.SF_pt_mn   = mn(SF_pt);                      R.SF_pt_sd   = sd(SF_pt);
R.SF_be_mn   = mn(SF_be);                      R.SF_be_sd   = sd(SF_be);
R.dSNR_pt_mn = mn(pct(SNR_pt,x_P));            R.dSNR_pt_sd = sd(pct(SNR_pt,x_P));
R.dSNR_be_mn = mn(pct(SNR_be,x_I));            R.dSNR_be_sd = sd(pct(SNR_be,x_I));
R.dSF_pt_mn  = mn(pct(SF_pt,1./(1+1./x_P)));   R.dSF_pt_sd  = sd(pct(SF_pt,1./(1+1./x_P)));
R.dSF_be_mn  = mn(pct(SF_be,1./(1+1./x_I)));   R.dSF_be_sd  = sd(pct(SF_be,1./(1+1./x_I)));
R.sat_pt_max = max(sat_pt,[],3);
R.sat_be_max = max(sat_be,[],3);

% ------------------------------------------------------------
% Console summary (rows: closed-form S/N)
% ------------------------------------------------------------

print_table('Point source, pupil plane: Delta S/N (%) vs kappa_P*S/N (mean +/- std)', ...
    SNR_vec, names, R.dSNR_pt_mn, R.dSNR_pt_sd);
print_table('Point source, pupil plane: Delta S_F (%) vs 1/(1+1/(kappa_P*S/N)) (mean +/- std)', ...
    SNR_vec, names, R.dSF_pt_mn, R.dSF_pt_sd);
print_table('Bulls-eye, image plane: Delta S/N (%) vs kappa_I*S/N (mean +/- std)', ...
    SNR_vec, names, R.dSNR_be_mn, R.dSNR_be_sd);
print_table('Bulls-eye, image plane: Delta S_F (%) vs 1/(1+1/(kappa_I*S/N)) (mean +/- std)', ...
    SNR_vec, names, R.dSF_be_mn, R.dSF_be_sd);
print_table('Saturated hologram pixels (%), point source (worst realization)', ...
    SNR_vec, names, 100*R.sat_pt_max, []);
print_table('Saturated hologram pixels (%), bulls-eye (worst realization)', ...
    SNR_vec, names, 100*R.sat_be_max, []);

% ------------------------------------------------------------
% Figure
% ------------------------------------------------------------

figC = plot_comparison(R, names, figW, figH, fontSize);

if exportFigures
    exportgraphics(figC,'../outputs/comparison_2x2.pdf', ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(figC,'../outputs/comparison_2x2.png','Resolution',600);
end

if saveResults
    save('../outputs/comparison_results.mat', 'names','SNR_vec','nReals', ...
        'N','K','q_I','q_P','S','alpha','kappa_P','kappa_I','eta_all', ...
        'SF_pt','SNR_pt','SF_be','SNR_be','sat_pt','sat_be','R');
end

fprintf('\nDone.\n');

% ============================================================
% Local functions
% ============================================================

function g = offaxis_geometry(name, U_R, q, N)
% Off-axis geometry: one hologram, sideband at (+M/4,+M/4).
    g.name   = name;
    g.type   = 'offaxis';
    g.refs   = {U_R};
    g.M      = round(q*N);
    sb       = g.M/2 + 1 + g.M/4;
    g.idx    = (sb-N/2):(sb+N/2-1);
    [cx,cy]  = meshgrid(linspace(-(N-1)/2,(N-1)/2,N));
    g.circWin = hypot(cx,cy) <= (N-1)/2;
    g.Npad   = round(q*N);
end

function g = psrg_geometry(name, U_R, q, N)
% On-axis phase-shifting geometry: four holograms, piston shifts 0..3pi/2.
    g.name  = name;
    g.type  = 'psrg';
    d       = [0, pi/2, pi, 3*pi/2];
    g.refs  = arrayfun(@(dk) U_R*exp(-1j*dk), d, 'UniformOutput', false);
    g.Npad  = round(q*N);
end

function [dH0, sat] = record_registered(g, U_S, K, sigma_r, usePoisson, ell, enforce, g_AD, nBits)
% Record, register to the nominal grid, and remove the mean of frame 1.
    [~,~,dH] = record_holograms(U_S, g.refs, K, sigma_r, usePoisson, ell, ...
        enforce, g_AD, nBits);
    sat = mean(dH(:) == 2^nBits-1);
    dH  = register_pixel_centres(dH, K);
    f1  = dH(:,:,1);
    dH0 = dH - round(mean(f1(:)));
end

function snr = support_snr(e, n)
% Numerical S/N on a support from signal+noise samples e and
% reference-only samples n (same processing, same support).
    snr = (mean(abs(e).^2) - mean(abs(n).^2))/var(n);
end

function UPe = pupil_estimate(g, dH0, g_AD, pupilMask)
% Pupil-plane estimate for S_F (same path as in the illustrative scripts).
    N = size(dH0,1);
    switch g.name
        case 'IPRG'
            UPe = reconstruct_estimate_iprg(dH0/g_AD, g.M, g.idx, g.idx, pupilMask);
        case 'PPRG'
            UIe = reconstruct_estimate_pprg(dH0/g_AD, g.M, g.idx, g.idx, pupilMask);
            UPe = pad_to_pupil(UIe, g.Npad, pupilMask);
        case 'IP-PSRG'
            USe = reconstruct_estimate_psrg(dH0/g_AD, ones(N));
            UPe = pad_to_pupil(USe, g.Npad, pupilMask);
        case 'PP-PSRG'
            UPe = reconstruct_estimate_psrg(dH0/g_AD, pupilMask);
    end
end

function I = image_estimate(g, dH0, g_AD, pupilMask)
% Image-plane estimate (same path as Fig. 4(d) of the illustrative scripts).
    N = size(dH0,1);
    switch g.name
        case 'IPRG'      % pupil estimate -> zero-pad -> fft2c -> crop
            UPe = reconstruct_estimate_iprg(dH0/g_AD, g.M, g.idx, g.idx, pupilMask);
            I = pad_to_image(UPe, g.Npad);
        case 'PPRG'      % the sideband is the image
            I = reconstruct_estimate_pprg(dH0/g_AD, g.M, g.idx, g.idx, pupilMask);
        case 'IP-PSRG'   % the phase-shift estimate is the image
            I = reconstruct_estimate_psrg(dH0/g_AD, ones(N));
        case 'PP-PSRG'   % pupil estimate -> zero-pad -> fft2c -> crop
            UPe = reconstruct_estimate_psrg(dH0/g_AD, pupilMask);
            I = pad_to_image(UPe, g.Npad);
    end
end

function I = pad_to_image(UPe, Npad)
% Zero-pad the N x N pupil estimate to Npad, fft2c, crop N x N.
    N  = size(UPe,1);
    ip = centered_indices(Npad, N);
    Up = zeros(Npad);  Up(ip,ip) = UPe;
    t  = fft2c(Up);
    I  = t(ip,ip);
end

function UPe = pad_to_pupil(UIe, Npad, pupilMask)
% Zero-pad the N x N image-plane estimate to Npad, ifft2c, crop N x N, mask.
    N  = size(UIe,1);
    ip = centered_indices(Npad, N);
    Up = zeros(Npad);  Up(ip,ip) = UIe;
    t  = ifft2c(Up);
    UPe = pupilMask .* t(ip,ip);
end

function O = make_bullseye_target_phys(xO, yO, D)
    O = double(hypot(xO,yO) <= D/2) ...
      - double(hypot(xO,yO) <= 3*D/8) ...
      - double(hypot(xO,yO) <= D/4) ...
      + double(hypot(xO,yO) <= D/8);
end

function print_table(title, SNR_vec, names, mn, sd, nd)
% Console table: one row per SNR point, one column per geometry.
% nd = number of decimals (default 2).
    if nargin < 6, nd = 2; end
    fprintf('\n  %s\n', title);
    fprintf('  %8s', 'S/N');
    fprintf('  %16s', names{:});
    fprintf('\n');
    for si = 1:numel(SNR_vec)
        fprintf('  %8.3f', SNR_vec(si));
        for gi = 1:numel(names)
            if isempty(sd)
                fprintf('  %16.*f', nd+1, mn(gi,si));
            else
                fprintf('  %8.*f +/-%5.*f', nd, mn(gi,si), nd, sd(gi,si));
            end
        end
        fprintf('\n');
    end
end

function figC = plot_comparison(R, names, figW, figH, fontSize)
% 2x2 comparison figure, figW x figH inches, one shared legend on top.
% Top row: point source, pupil plane, vs kappa_P*S/N.
% Bottom row: bulls-eye, image plane, vs kappa_I*S/N.
    clr = [0.00 0.45 0.74;    % IPRG     blue
           0.85 0.33 0.10;    % PPRG     red
           0.47 0.67 0.19;    % IP-PSRG  green
           0.00 0.00 0.00];   % PP-PSRG  black
    mk  = {'o','s','d','^'};
    st.msz = 4;  st.lw = 1.0;  st.lweb = 0.6;  st.cs = 2;
    st.gry = [0.55 0.55 0.55];  st.clr = clr;  st.mk = mk;
    xlP = [10^floor(log10(min(R.x_P(:)))), 10^ceil(log10(max(R.x_P(:))))];
    xlI = [10^floor(log10(min(R.x_I(:)))), 10^ceil(log10(max(R.x_I(:))))];
    labP = '$\kappa_{\mathrm{P}}\,S/N$';
    labI = '$\kappa_{\mathrm{I}}\,S/N$';

    figC = make_figure(figW, figH);
    tC   = tiledlayout(figC,2,2,'TileSpacing','compact','Padding','compact');

    % (a) S_F, point source, pupil plane
    ax = nexttile(tC);
    [h, hTh] = sf_panel(ax, R.x_P, R.SF_pt_mn, R.SF_pt_sd, xlP, st);
    xlabel(ax,labP);  ylabel(ax,'$S_{\mathrm{F}}$');  panel_label(ax,'(a)');

    % (b) Delta S/N, point source, pupil plane
    ax = nexttile(tC);
    pct_panel(ax, R.x_P, R.dSNR_pt_mn, R.dSNR_pt_sd, xlP, st);
    xlabel(ax,labP);  ylabel(ax,'$\Delta S/N$ (\%)');  panel_label(ax,'(b)');

    % (c) S_F, extended object, image plane
    ax = nexttile(tC);
    sf_panel(ax, R.x_I, R.SF_be_mn, R.SF_be_sd, xlI, st);
    xlabel(ax,labI);  ylabel(ax,'$S_{\mathrm{F}}$');  panel_label(ax,'(c)');

    % (d) Delta S/N, extended object, image plane
    ax = nexttile(tC);
    pct_panel(ax, R.x_I, R.dSNR_be_mn, R.dSNR_be_sd, xlI, st);
    xlabel(ax,labI);  ylabel(ax,'$\Delta S/N$ (\%)');  panel_label(ax,'(d)');

    % One legend for all panels, above the tiles
    lgd = legend([h, hTh], [names, {'Analytical'}], 'Orientation','horizontal', ...
        'NumColumns', numel(names)+1, 'FontSize', fontSize-1, 'Box','off');
    lgd.Layout.Tile = 'north';
end

function [h, hTh] = sf_panel(ax, x, mn, sd, xl, st)
% S_F panel: analytical 1/(1+1/x) and numerical mean +/- std.
    xf  = logspace(log10(xl(1)), log10(xl(2)), 300);
    hTh = semilogx(ax, xf, 1./(1+1./xf), '-', 'Color', st.gry, 'LineWidth', st.lw);
    hold(ax,'on');
    h = plot_markers(ax, x, mn, sd, st, 'none');
    style_axes(ax, xl);
    ylim(ax,[0 1]);  yticks(ax,0:0.2:1);
end

function pct_panel(ax, x, mn, sd, xl, st)
% Percentage-error panel with a zero reference line.
    yline(ax, 0, '--', 'Color', st.gry, 'LineWidth', 0.5);
    hold(ax,'on');
    plot_markers(ax, x, mn, sd, st, '-');
    style_axes(ax, xl);
end

function h = plot_markers(ax, x, mn, sd, st, ls)
% Mean +/- std markers, one series per geometry.
    h = gobjects(1,size(x,1));
    for gi = 1:size(x,1)
        h(gi) = errorbar(ax, x(gi,:), mn(gi,:), sd(gi,:), ...
            'LineStyle', ls, 'Marker', st.mk{gi}, 'Color', st.clr(gi,:), ...
            'MarkerFaceColor', 'w', 'MarkerSize', st.msz, ...
            'LineWidth', st.lweb, 'CapSize', st.cs);
    end
end

function style_axes(ax, xl)
% Common log-x axis styling.
    set(ax, 'XScale','log', 'XMinorTick','on', 'YMinorTick','on', ...
        'TickDir','in', 'Layer','top');
    xlim(ax, xl);  xticks(ax, 10.^(log10(xl(1)):log10(xl(2))));
    grid(ax,'on');  ax.MinorGridLineStyle = ':';  box(ax,'on');
end

function panel_label(ax, str)
% Panel label above the top-left corner (same as the illustrative figures).
    text(ax, 0, 1.02, str, 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end
