% ============================================================
% script_PPPSRG_illustrative.m
%
% Illustrative examples for the on-axis PP-PSRG (Section 6.5).
% (Pupil-Plane Phase-Shifting Recording Geometry)
%
% Produces four figures:
%   Figure 1 (3x2): point-source holograms and pupil estimate
%   Figure 2 (2x2): point-source truth vs. estimate, pupil plane
%   Figure 3 (3x2): bulls-eye holograms and pupil estimate
%   Figure 4 (2x2): bulls-eye reflectance and images
%
% KEY CHARACTERISTICS:
%   - Pupil P generated on fine Nfine×Nfine grid at spacing delta_f = delta/K.
%   - FPA in PUPIL plane; U_S_fine = P_fine (no propagation needed).
%   - On-axis uniform reference on fine grid; four piston-shifted measurements.
%   - pixel_bin on hologram irradiance -> N×N detector holograms.
%   - Downsampling of P_fine (every K-th sample) -> N×N pupil truth.
%   - S_F computed in PUPIL plane between P_coarse and USe.
%   - Fig 4 images via fft2c of pupil fields (Fourier = image plane).
%   - Reference-only noise estimate from a fifth recording (S=0).
%
% All functions are in ../functions/.
% ============================================================

clear; close all; clc;
rng(1);

addpath('../functions');

% ============================================================
% Parameters
% ============================================================
N      = 256;         % number of detector pixels per side
K      = 2;           % pixel binning factor
D      = 0.30;        % exit-pupil diameter [m]
delta  = D/N;         % detector pixel spacing [m]
delta_f= delta/K;     % sub-pixel (fine) spacing [m]
Nfine  = K*N;         % fine-grid size
wvl    = 1.0e-6;      % wavelength [m]
kwave  = 2*pi/wvl;    % angular wavenumber [rad/m]
q_P    = 3;           % pixels per speckle at pupil plane
z_P    = q_P*D^2/(N*wvl); % speckle = q_P*delta, same formula as z_I [m]
S      = 1;           % amplitude splits
sigma_r= 100;         % read-noise std [pe]
ell    = 100e3;        % well depth [pe]
m_R    = 0.75*ell;    % reference mean pe
m_S_pt = 1;           % point-source signal mean pe
m_S_be = 4;           % bulls-eye signal mean pe
nBits  = 14;          % A/D bits
g_AD   = (2^nBits-1)/ell;
phi_rms   = 1.0;
phaseCorr = D/10;

sigma_n2  = m_R + sigma_r^2 + (1/12)/g_AD^2;
sigma_PS2 = sigma_n2/4;

print_run_header('ON-AXIS PP-PSRG ILLUSTRATIVE EXAMPLES',N,K,'q_P',q_P,D,wvl,'z_P',z_P,delta,nBits,g_AD,sigma_n2);

% ============================================================
% Fine pupil-plane grid (Nfine×Nfine at delta_f)
% ============================================================
x1f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[x1f,y1f] = meshgrid(x1f_vec, x1f_vec);
r1f = hypot(x1f,y1f);
pupilMask_fine = double(r1f <= D/2);

phi_fine = make_isoplanatic_phase(Nfine, delta_f, phaseCorr, phi_rms, pupilMask_fine);
P_fine   = pupilMask_fine .* exp(1j*phi_fine);

% Coarse pupil (N×N): alias-free decimation
P_coarse   = P_fine(1:K:end, 1:K:end);
pupilMask  = pupilMask_fine(1:K:end, 1:K:end);
x1_vec     = x1f_vec(1:K:end);

% ============================================================
% Uniform reference on fine grid (four piston shifts)
% ============================================================
A_R = sqrt(m_R);
delta_piston = [0, pi/2, pi, 3*pi/2];
U_R_cell = cell(1,4);
for kk = 1:4
    U_R_cell{kk} = A_R*exp(-1j*delta_piston(kk))*ones(Nfine);
end

% Closed-form SNR per FPA pixel (m_S = full-FPA mean), with pixel-aperture
% efficiency eta, consistent with the off-axis geometries.
snr_theory = @(mS,eta) eta*(4/S)*m_R*mS/sigma_n2;
% S_F is evaluated on the pupil support, which holds all of the signal but
% only pi/4 of the FPA pixels. SNR on the pupil support = kappa * SNR.
kappa = 4/pi;

% ============================================================
% Reference-only noise recording
% ============================================================
U_ref_only = reference_only_pppsrg(U_R_cell, sigma_r, ell, g_AD, nBits, K);

% ============================================================
% Figure formatting
% ============================================================
figW=6.5; figH3=6.5; figH2=5.0; fontSize=10;

set(groot,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultColorbarTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',fontSize,'defaultTextFontSize',fontSize, ...
    'defaultColorbarFontSize',fontSize,'defaultAxesLineWidth',0.5);

uPup      = x1_vec/D;
uPix      = (0:N-1)-N/2;
uPix_fine = (0:Nfine-1)-Nfine/2;

labPx='$x_1/D$'; labPy='$y_1/D$';
labNx='$x_1/p$'; labNy='$y_1/p$';
cbAmpUP  ='$|U_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUP  ='$\angle U_\mathrm{P}\;(\mathrm{rad})$';
cbAmpUSe ='$|\hat{U}_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUSe ='$\angle\hat{U}_\mathrm{P}\;(\mathrm{rad})$';
cbHolo0  ='$d_\mathrm{H}^{(\delta=0)}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloPi2='$d_\mathrm{H}^{(\delta=\pi/2)}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloN0 ='$d_\mathrm{H}^{+(\delta=0)}-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbHoloNPi2='$d_\mathrm{H}^{+(\delta=\pi/2)}-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbAmpRO  ='$|R_\mathrm{O}|\;(\mathrm{norm.})$';
cbPhzRO  ='$\angle R_\mathrm{O}\;(\mathrm{rad})$';
cbIrrUI  ='$|U_\mathrm{I}|^2\;(\mathrm{norm.})$';
cbIrrUIe ='$|\hat{U}_\mathrm{I}|^2\;(\mathrm{norm.})$';

% ============================================================
% PART I: POINT-SOURCE OBJECT
% ============================================================

% Signal = pupil field directly (FPA in pupil plane)
U_S_pt_fine = sqrt(m_S_pt)*P_fine/sqrt(mean(abs(P_fine(:)).^2));

[dH_pt_ideal, dH_pt] = record_holograms_pppsrg(U_S_pt_fine, U_R_cell, S, ...
    sigma_r, ell, g_AD, nBits, K);
holoMaxDN_pt = max(dH_pt(:));
% Binned values sit at pixel centres; register them to the nominal grid
dH_pt_ideal = register_pixel_centres(dH_pt_ideal, K);
dH_pt       = register_pixel_centres(dH_pt, K);

dH0_pt_ideal = dH_pt_ideal - round(mean(dH_pt_ideal(:,:,1),'all'));
dH0_pt       = dH_pt       - round(mean(dH_pt(:,:,1),'all'));

USe_pt0 = reconstruct_estimate_psrg(dH0_pt_ideal/g_AD, pupilMask);
USe_pt  = reconstruct_estimate_psrg(dH0_pt/g_AD,       pupilMask);

SF_pt0 = field_strehl(P_coarse, USe_pt0, pupilMask);
SF_pt  = field_strehl(P_coarse, USe_pt,  pupilMask);
eta_pt = pixel_eta(U_S_pt_fine.*conj(U_R_cell{1}), K);
SNR_pt = snr_theory(m_S_pt, eta_pt);

% Noise stats: signal+noise and reference-only estimates must cover the
% SAME region. U_ref_only is unmasked, so use an unmasked estimate here.
USe_pt_full = reconstruct_estimate_psrg(dH0_pt/g_AD, ones(N));
nz_pt = pppsrg_noise_stats(USe_pt_full, U_ref_only, pupilMask);
psrg_report('Point source', U_S_pt_fine, U_R_cell, S, K, m_S_pt, m_R, holoMaxDN_pt, ...
    nBits, eta_pt, nz_pt, SNR_pt, sigma_PS2, kappa, SF_pt0, SF_pt);

USe_pt_a = align_global_phase(P_coarse, USe_pt, pupilMask);

fig1 = make_figure(figW,figH3);
t1   = tiledlayout(fig1,3,2,'TileSpacing','compact','Padding','compact');
cH_pt = max(abs(dH0_pt_ideal(:,:,1:2)),[],'all');
show_holo(t1,uPix,uPix,dH0_pt_ideal(:,:,1),cH_pt,labNx,labNy,cbHolo0,  '(a)');
show_holo(t1,uPix,uPix,dH0_pt_ideal(:,:,2),cH_pt,labNx,labNy,cbHoloPi2,'(b)');
show_holo(t1,uPix,uPix,dH0_pt(:,:,1),      cH_pt,labNx,labNy,cbHoloN0, '(c)');
show_holo(t1,uPix,uPix,dH0_pt(:,:,2),      cH_pt,labNx,labNy,cbHoloNPi2,'(d)');
show_amp  (t1,uPix,uPix,abs(USe_pt_a),pupilMask,labNx,labNy,cbAmpUSe,'(e)');
show_phase(t1,uPix,uPix,USe_pt_a,     pupilMask,labNx,labNy,cbPhzUSe,'(f)');

fig2 = make_figure(figW,figH2);
t2   = tiledlayout(fig2,2,2,'TileSpacing','compact','Padding','compact');
show_amp  (t2,uPup,uPup,abs(P_coarse), pupilMask,labPx,labPy,cbAmpUP, '(a)');
show_phase(t2,uPup,uPup,P_coarse,      pupilMask,labPx,labPy,cbPhzUP, '(b)');
show_amp  (t2,uPup,uPup,abs(USe_pt_a), pupilMask,labPx,labPy,cbAmpUSe,'(c)');
show_phase(t2,uPup,uPup,USe_pt_a,      pupilMask,labPx,labPy,cbPhzUSe,'(d)');

% ============================================================
% PART II: EXTENDED BULLS-EYE OBJECT
% ============================================================

% Object and reflectance on fine grid
[xOf,yOf] = meshgrid((-Nfine/2:Nfine/2-1)*delta_f);
O_fine   = make_bullseye_target_phys(xOf, yOf, D);
phz_fine = -pi + 2*pi*rand(Nfine);
R_O_fine = O_fine.*exp(1j*phz_fine);

% Propagate object -> pupil, apply aberration
U_Pm_fine = fresnel_prop_object_to_pupil(R_O_fine, D, delta_f, wvl, z_P);
U_P_be_fine = P_fine .* U_Pm_fine;
U_S_be_fine = sqrt(m_S_be)*U_P_be_fine/sqrt(mean(abs(U_P_be_fine(:)).^2));

% Coarse for S_F truth
U_P_be_coarse = U_P_be_fine(1:K:end, 1:K:end);

[dH_be_ideal, dH_be] = record_holograms_pppsrg(U_S_be_fine, U_R_cell, S, ...
    sigma_r, ell, g_AD, nBits, K);
holoMaxDN_be = max(dH_be(:));
dH_be_ideal = register_pixel_centres(dH_be_ideal, K);
dH_be       = register_pixel_centres(dH_be, K);

dH0_be_ideal = dH_be_ideal - round(mean(dH_be_ideal(:,:,1),'all'));
dH0_be       = dH_be       - round(mean(dH_be(:,:,1),'all'));

USe_be0 = reconstruct_estimate_psrg(dH0_be_ideal/g_AD, pupilMask);
USe_be  = reconstruct_estimate_psrg(dH0_be/g_AD,       pupilMask);

SF_be0 = field_strehl(U_P_be_coarse, USe_be0, pupilMask);
SF_be  = field_strehl(U_P_be_coarse, USe_be,  pupilMask);
eta_be = pixel_eta(U_S_be_fine.*conj(U_R_cell{1}), K);
SNR_be = snr_theory(m_S_be, eta_be);

USe_be_full = reconstruct_estimate_psrg(dH0_be/g_AD, ones(N));
nz_be = pppsrg_noise_stats(USe_be_full, U_ref_only, pupilMask);
psrg_report('Bulls-eye', U_S_be_fine, U_R_cell, S, K, m_S_be, m_R, holoMaxDN_be, ...
    nBits, eta_be, nz_be, SNR_be, sigma_PS2, kappa, SF_be0, SF_be);

USe_be_a = align_global_phase(U_P_be_coarse, USe_be, pupilMask);

fig3 = make_figure(figW,figH3);
t3   = tiledlayout(fig3,3,2,'TileSpacing','compact','Padding','compact');
cH_be = max(abs(dH0_be_ideal(:,:,1:2)),[],'all');
show_holo(t3,uPix,uPix,dH0_be_ideal(:,:,1),cH_be,labNx,labNy,cbHolo0,  '(a)');
show_holo(t3,uPix,uPix,dH0_be_ideal(:,:,2),cH_be,labNx,labNy,cbHoloPi2,'(b)');
show_holo(t3,uPix,uPix,dH0_be(:,:,1),      cH_be,labNx,labNy,cbHoloN0, '(c)');
show_holo(t3,uPix,uPix,dH0_be(:,:,2),      cH_be,labNx,labNy,cbHoloNPi2,'(d)');
show_amp  (t3,uPix,uPix,abs(USe_be_a),pupilMask,labNx,labNy,cbAmpUSe,'(e)');
show_phase(t3,uPix,uPix,USe_be_a,     pupilMask,labNx,labNy,cbPhzUSe,'(f)');

% Fig 4: reflectance (a)(b), truth image (c), estimated image (d)
% Truth: Fresnel propagate fine-grid pupil -> fine image (Nfine×Nfine)
U_I_be_fine = fresnel_prop_pupil_to_image(U_P_be_fine, D, delta_f, wvl, z_P);
Itru = abs(U_I_be_fine).^2;  Itru = Itru/max(Itru(:));
% Estimate: zero-pad to q_P*N, fft2c, crop central N×N
Npad_est = q_P*N;
USe_be_pad = zeros(Npad_est);
ip_est = centered_indices(Npad_est, N);
USe_be_pad(ip_est,ip_est) = USe_be_a;
Iest_pad = abs(fft2c(USe_be_pad)).^2;
ic   = centered_indices(Npad_est, N);
Iest = Iest_pad(ic,ic); Iest = Iest/max(Iest(:));
objMask_fine = double(O_fine ~= 0);

fig4 = make_figure(figW,figH2);
t4   = tiledlayout(fig4,2,2,'TileSpacing','compact','Padding','compact');
labWx='$x/W$'; labWy='$y/W$';
uObj_fine = x1f_vec/D;   % fine grid normalized by W=D
uObj      = x1_vec/D;    % coarse grid normalized by W=D
show_amp  (t4,uObj_fine,uObj_fine,abs(R_O_fine),[],        labWx,labWy,cbAmpRO,'(a)');
show_phase(t4,uObj_fine,uObj_fine,R_O_fine,    objMask_fine,labWx,labWy,cbPhzRO,'(b)');
show_amp  (t4,uObj_fine,uObj_fine,Itru,[],labWx,labWy,cbIrrUI, '(c)');
show_amp  (t4,uObj,     uObj,     Iest,[],labWx,labWy,cbIrrUIe,'(d)');

exportgraphics(fig1,'../outputs/PPPSRG_pt_holograms.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig1,'../outputs/PPPSRG_pt_holograms.png','Resolution',600);
exportgraphics(fig2,'../outputs/PPPSRG_pt_comparison.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig2,'../outputs/PPPSRG_pt_comparison.png','Resolution',600);
exportgraphics(fig3,'../outputs/PPPSRG_bulls_holograms.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig3,'../outputs/PPPSRG_bulls_holograms.png','Resolution',600);
exportgraphics(fig4,'../outputs/PPPSRG_bulls_comparison.pdf', ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig4,'../outputs/PPPSRG_bulls_comparison.png','Resolution',600);

% ============================================================
% Local functions
% ============================================================

function [dH_ideal, dH] = record_holograms_pppsrg(U_S_fine, U_R_cell, S, ...
    sigma_r, ell, g_AD, nBits, K)
    Nfine4 = size(U_S_fine,1);
    Nd     = Nfine4/K;
    dH_ideal = zeros(Nd,Nd,4);
    dH       = zeros(Nd,Nd,4);
    usePoisson = license('test','statistics_toolbox');
    for kk = 1:4
        mH_fine = abs(U_S_fine/S + U_R_cell{kk}).^2;
        mH = pixel_bin(mH_fine, K);
        if usePoisson
            mHn = double(poissrnd(mH)) + sigma_r*randn(Nd);
        else
            mHn = mH + sqrt(mH).*randn(Nd) + sigma_r*randn(Nd);
        end
        mHn = min(mHn, ell);
        dH_ideal(:,:,kk) = round(g_AD*mH);
        dH(:,:,kk)       = round(g_AD*mHn);
    end
    dH_ideal = min(dH_ideal, 2^nBits-1);
    dH       = min(dH,       2^nBits-1);
end

function U_N = reference_only_pppsrg(U_R_cell, sigma_r, ell, g_AD, nBits, K)
    Nfine4   = size(U_R_cell{1},1);
    U_S_zero = zeros(Nfine4);
    [~, dH_n] = record_holograms_pppsrg(U_S_zero, U_R_cell, 1, ...
        sigma_r, ell, g_AD, nBits, K);
    dH_n  = register_pixel_centres(dH_n, K);
    Nd    = size(dH_n,1);
    dH0_n = dH_n - round(mean(dH_n(:,:,1),'all'));
    U_N   = reconstruct_estimate_psrg(dH0_n/g_AD, ones(Nd));
end

function nz = pppsrg_noise_stats(USe_sn, U_N, mask)
    % SNR computed over full grid (consistent with 2018 paper)
    s = USe_sn(:);  n = U_N(:);
    nz.powerSN     = mean(abs(s).^2);
    nz.powerN      = mean(abs(n).^2);
    nz.varNoise    = var(n);
    nz.signalPower = nz.powerSN - nz.powerN;
    nz.SNR         = nz.signalPower / nz.varNoise;
end

function O = make_bullseye_target_phys(xO, yO, D)
    O = double(hypot(xO,yO) <= D/2) ...
      - double(hypot(xO,yO) <= 3*D/8) ...
      - double(hypot(xO,yO) <= D/4) ...
      + double(hypot(xO,yO) <= D/8);
end

function psrg_report(name,U_S,U_R_cell,S,K,mS,mR,holoMaxDN,nBits,eta,nz,snrClosed,sigma_PS2,kappa,SF0,SF)
% Standard console block (PRINT_OBJECT_REPORT) for the phase-shifting geometries.
% Noise and signal are per pixel of the signal estimate [pe^2].
    r.mS = mS;  r.mS_meas = mean(pixel_bin(abs(U_S).^2,K),'all');
    r.mR = mR;  r.mR_meas = mean(pixel_bin(abs(U_R_cell{1}).^2,K),'all');
    mx = 0;
    for kk = 1:numel(U_R_cell)
        mH = pixel_bin(abs(U_S/S + U_R_cell{kk}).^2,K);
        mx = max(mx, max(mH(:)));
    end
    r.holoMaxPe = mx;  r.holoMaxDN = holoMaxDN;  r.maxDN = 2^nBits-1;
    r.eta = eta;
    r.noisePred = sigma_PS2;         r.noiseMeas = nz.varNoise;
    r.nEff      = numel(U_S)/K^2;    % noise stats over the full N x N grid
    r.sigPred   = eta*mR*mS/S;       r.sigMeas   = nz.signalPower;
    r.snrClosed = snrClosed;  r.snrNum = nz.SNR;
    r.kappa = kappa;
    r.SF0 = SF0;  r.SF = SF;
    print_object_report(name,r);
end
