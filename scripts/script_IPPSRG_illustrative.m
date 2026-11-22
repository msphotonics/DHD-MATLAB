% ============================================================
% script_IPPSRG_illustrative.m
%
% Illustrative examples for the on-axis IP-PSRG.
% (Image-Plane Phase-Shifting Recording Geometry)
%
% Produces two figures (3x2), saved as IPPSRG_pt and IPPSRG_bulls:
%   Figure 1: point-source (a),(b) holograms for delta = 0, (c),(d) holograms for delta = pi,
%             (e),(f) truth vs. estimated pupil-plane phase
%   Figure 3: bulls-eye (a),(b) holograms for delta = 0, (c),(d) holograms for delta = pi,
%             (e),(f) truth vs. estimated image-plane irradiance
%
% KEY CHARACTERISTICS:
%   - Pupil P generated on fine Nfine×Nfine grid at spacing delta_f = delta/K.
%   - Fresnel propagation on fine grid -> U_I_fine (Nfine×Nfine).
%   - pixel_bin on hologram irradiance -> N×N detector holograms.
%   - Downsampling of P_fine (every K-th sample) -> N×N pupil truth.
%   - Reconstruction: reconstruct_estimate_psrg -> USe (N×N image plane).
%   - ifft2c(USe) -> U_P_est (N×N pupil plane).
%   - S_F computed in PUPIL plane between P_coarse and U_P_est.
%   - Reference-only noise estimate from a fifth recording (U_S = 0).
%
% All functions are in ../functions/.
% ============================================================

clear; close all; clc;
rng(1);

% Run from this script's folder, so ../functions and ../outputs resolve
% wherever MATLAB's current folder happens to be
cd(fileparts(mfilename('fullpath')));
addpath('../functions');

% ============================================================
% Parameters
% ============================================================
N      = 256;         % number of detector pixels per side
K      = 2;           % pixel binning factor (sub-pixels per detector pixel)
D      = 0.30;        % exit-pupil diameter [m]
delta  = D/N;         % detector pixel spacing [m]
delta_f= delta/K;     % sub-pixel (fine) spacing [m]
Nfine  = K*N;         % fine-grid size
wvl    = 1.0e-6;      % wavelength [m]
kwave  = 2*pi/wvl;    % angular wavenumber [rad/m]
q_I    = 3;           % IP-PSRG quotient (detector pixels per PSF)
ifov   = wvl/D/q_I;   % pixel field of view [rad]
z_I    = D/(N*ifov);  % propagation distance [m]
B      = 1;           % amplitude splits (B in the chapter; 1 = sequential)
sigma_r= 100;         % read-noise std [pe]
ell    = 100e3;        % well depth [pe]
m_R    = 0.75*ell;    % reference mean pe
m_S_pt = 1;           % point-source signal mean pe
m_S_be = 4;           % bulls-eye signal mean pe
nBits  = 14;          % A/D bits
g_AD   = (2^nBits-1)/ell; % A/D gain [DN/pe]
phi_rms   = 1.0;      % phase aberration rms [rad]
phaseCorr = D/10;     % correlation length [m]

sigma_n2  = m_R + sigma_r^2 + (1/12)/g_AD^2;
sigma_PS2 = sigma_n2/4;

print_run_header('ON-AXIS IP-PSRG ILLUSTRATIVE EXAMPLES',N,K,'q_I',q_I,D,wvl,'z_I',z_I,delta,nBits,g_AD,sigma_n2);

% ============================================================
% Fine pupil-plane grid (Nfine×Nfine at delta_f)
% ============================================================
x1f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[x1f,y1f] = meshgrid(x1f_vec, x1f_vec);
r1f = hypot(x1f,y1f);
pupilMask_fine = double(r1f <= D/2);

% Isoplanatic phase on fine grid
phi_fine = make_isoplanatic_phase(Nfine, delta_f, phaseCorr, phi_rms, pupilMask_fine);
P_fine   = pupilMask_fine .* exp(1j*phi_fine);

% Coarse pupil (N×N): downsample P_fine by taking every K-th sample
% Valid because pupil is bandlimited to D — alias-free decimation
P_coarse    = P_fine(1:K:end, 1:K:end);
pupilMask   = pupilMask_fine(1:K:end, 1:K:end);
x1_vec      = x1f_vec(1:K:end);
[x1,y1]    = meshgrid(x1_vec, x1_vec);

% ============================================================
% Fine image-plane grid (Nfine×Nfine at delta_f)
% ============================================================
x2f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[x2f,y2f] = meshgrid(x2f_vec, x2f_vec);

% Coarse image-plane grid (N×N at delta)
x2_vec = x2f_vec(1:K:end);
[x2,y2] = meshgrid(x2_vec, x2_vec);

% ============================================================
% Reference wave on fine grid: positive quadratic phase (Eq. 28)
% ============================================================
U_R_base = (1/1j)*exp(1j*kwave*z_I) ...
           .*exp(+1j*kwave/(2*z_I)*(x2f.^2+y2f.^2));
U_R_base = sqrt(m_R)*U_R_base/sqrt(mean(abs(U_R_base(:)).^2));

delta_piston = [0, pi/2, pi, 3*pi/2];
U_R_cell = cell(1,4);
for kk = 1:4
    U_R_cell{kk} = U_R_base*exp(-1j*delta_piston(kk));
end

% Closed-form SNR per image-plane pixel, with pixel-aperture efficiency eta
snr_theory = @(mS,eta) eta*(4/B)*m_R*mS/sigma_n2;
% S_F is evaluated in the PUPIL plane after the pupil mask, which keeps only
% the pupil band (a fraction pi/(4 q_I^2) of the image-plane noise) and all
% of the signal. SNR on the pupil support = kappa * SNR.
kappa = 4*q_I^2/pi;

% Truth phase alignment: remove quadratic phase to match USe
quad_phase_remove = @(U_S) exp(1j*pi/2)*U_S ...
    .*exp(-1j*kwave/(2*z_I)*(x2.^2+y2.^2));

% ============================================================
% Reference-only noise recording
% ============================================================
U_ref_only = reference_only_ippsrg(U_R_cell, sigma_r, ell, g_AD, nBits, K);

% ============================================================
% Figure formatting
% ============================================================
figW=6.5; figH3=6.5; figH2=5.0; figH1=2.7; fontSize=10;

set(groot,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultColorbarTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',fontSize,'defaultTextFontSize',fontSize, ...
    'defaultColorbarFontSize',fontSize,'defaultAxesLineWidth',0.5);

uPix = (0:N-1)-N/2;
uPup = x1_vec/D;

labNx='$x_2/p$'; labNy='$y_2/p$';
labPx='$x_1/D$'; labPy='$y_1/D$';
cbHolo0  ='$d_\mathrm{H}^{(\delta=0)}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloPi2='$d_\mathrm{H}^{(\delta=\pi/2)}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloN0 ='$d_\mathrm{H}^{+(\delta=0)}-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbHoloNPi2='$d_\mathrm{H}^{+(\delta=\pi/2)}-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbHoloPi ='$d_\mathrm{H}^{(\delta=\pi)}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloNPi='$d_\mathrm{H}^{+(\delta=\pi)}-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbAmpUSe ='$|\hat{U}_\mathrm{S}|\;(\mathrm{norm.})$';
cbPhzUSe ='$\angle\hat{U}_\mathrm{S}\;(\mathrm{rad})$';
cbAmpUP  ='$|U_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUP  ='$\angle U_\mathrm{P}\;(\mathrm{rad})$';
cbAmpUPe ='$|\hat{U}_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUPe ='$\angle\hat{U}_\mathrm{P}\;(\mathrm{rad})$';
cbAmpRO  ='$|R_\mathrm{O}|\;(\mathrm{norm.})$';
cbPhzRO  ='$\angle R_\mathrm{O}\;(\mathrm{rad})$';
cbIrrUS  ='$|U_\mathrm{S}|^2\;(\mathrm{norm.})$';
cbIrrUSe ='$|\hat{U}_\mathrm{S}|^2\;(\mathrm{norm.})$';

% ============================================================
% PART I: POINT-SOURCE OBJECT
% ============================================================

% Propagate fine pupil -> fine image plane
U_I_pt_fine = fresnel_prop_pupil_to_image(P_fine, D, delta_f, wvl, z_I);
U_S_pt_fine = sqrt(m_S_pt)*U_I_pt_fine/sqrt(mean(abs(U_I_pt_fine(:)).^2));

% Downsample for display and S_F
U_S_pt = U_S_pt_fine(1:K:end, 1:K:end);

[dH_pt_ideal, dH_pt] = record_holograms_ippsrg(U_S_pt_fine, U_R_cell, B, ...
    sigma_r, ell, g_AD, nBits, K);
holoMaxDN_pt = max(dH_pt(:));
% Binned values sit at pixel centres; register them to the nominal grid
dH_pt_ideal = register_pixel_centres(dH_pt_ideal, K);
dH_pt       = register_pixel_centres(dH_pt, K);

dH0_pt_ideal = dH_pt_ideal - round(mean(dH_pt_ideal(:,:,1),'all'));
dH0_pt       = dH_pt       - round(mean(dH_pt(:,:,1),'all'));

% Reconstruct image-plane estimate (N×N)
USe_pt0 = reconstruct_estimate_psrg(dH0_pt_ideal/g_AD, ones(N));
USe_pt  = reconstruct_estimate_psrg(dH0_pt/g_AD,       ones(N));

% Pupil estimate: zero-pad USe to q_I*N before ifft2c to correct magnification
% USe is N×N at spacing delta_I = q_I*delta; ifft2c of q_I*N grid gives spacing delta
Npad = q_I*N;
ip   = centered_indices(Npad, N);
USe_pt0_pad = zeros(Npad); USe_pt0_pad(ip,ip) = USe_pt0;
USe_pt_pad  = zeros(Npad); USe_pt_pad(ip,ip)  = USe_pt;
tmp0 = ifft2c(USe_pt0_pad); U_P_est_pt0 = pupilMask .* tmp0(ip,ip);
tmp  = ifft2c(USe_pt_pad);  U_P_est_pt  = pupilMask .* tmp(ip,ip);

% S_F in pupil plane: truth = P_coarse, estimate = U_P_est
SF_pt0 = field_strehl(P_coarse, U_P_est_pt0, pupilMask);
SF_pt  = field_strehl(P_coarse, U_P_est_pt,  pupilMask);
% eta acts on the cross term U_S*conj(U_R): the spherical reference removes
% the image-plane quadratic phase before pixel binning.
eta_pt = pixel_eta(U_S_pt_fine.*conj(U_R_base), K);
SNR_pt = snr_theory(m_S_pt, eta_pt);

nz_pt = ippsrg_noise_stats(USe_pt, U_ref_only);
psrg_report('Point source', U_S_pt_fine, U_R_cell, B, K, m_S_pt, m_R, holoMaxDN_pt, ...
    nBits, eta_pt, nz_pt, SNR_pt, sigma_PS2, kappa, SF_pt0, SF_pt);

% Fig 1 (a-d): delta = 0 and pi holograms
fig1 = make_figure(figW,figH3);
t1   = tiledlayout(fig1,3,2,'TileSpacing','compact','Padding','compact');
cH_pt = max(abs(dH0_pt_ideal(:,:,[1 3])),[],'all');   % delta = 0 and pi
ax_ = show_holo(t1,uPix,uPix,dH0_pt_ideal(:,:,1),cH_pt,labNx,labNy,cbHolo0, '(a)');
add_zoom_inset(ax_,uPix,uPix,dH0_pt_ideal(:,:,1),[-8 8],[-8 8],[40 124],[40 124]);
ax_ = show_holo(t1,uPix,uPix,dH0_pt(:,:,1),      cH_pt,labNx,labNy,cbHoloN0,'(b)');
add_zoom_inset(ax_,uPix,uPix,dH0_pt(:,:,1),[-8 8],[-8 8],[40 124],[40 124]);
USe_pt_disp = align_global_phase(quad_phase_remove(U_S_pt), USe_pt, true(N));
ax_ = show_holo(t1,uPix,uPix,dH0_pt_ideal(:,:,3),cH_pt,labNx,labNy,cbHoloPi, '(c)');
add_zoom_inset(ax_,uPix,uPix,dH0_pt_ideal(:,:,3),[-8 8],[-8 8],[40 124],[40 124]);
ax_ = show_holo(t1,uPix,uPix,dH0_pt(:,:,3),      cH_pt,labNx,labNy,cbHoloNPi,'(d)');
add_zoom_inset(ax_,uPix,uPix,dH0_pt(:,:,3),[-8 8],[-8 8],[40 124],[40 124]);

% Fig 1 (e,f): pupil truth vs pupil estimate
U_P_est_pt_a = align_global_phase(P_coarse, U_P_est_pt, pupilMask);
show_phase(t1,uPup,uPup,P_coarse,         pupilMask,labPx,labPy,cbPhzUP, '(e)');
show_phase(t1,uPup,uPup,U_P_est_pt_a,     pupilMask,labPx,labPy,cbPhzUPe,'(f)');

% ============================================================
% PART II: EXTENDED BULLS-EYE OBJECT
% ============================================================

% Object on fine grid
[xOf,yOf] = meshgrid((-Nfine/2:Nfine/2-1)*delta_f);
O_fine  = make_bullseye_target_phys(xOf, yOf, D);
phz_fine = -pi + 2*pi*rand(Nfine);
R_O_fine = O_fine.*exp(1j*phz_fine);

% Physical propagation: object -> pupil (collimate) -> image (focus)
% z_P = z_I for unit magnification
U_Pm_be_fine = fresnel_prop_object_to_pupil(R_O_fine, D, delta_f, wvl, z_I);
U_P_be_fine  = P_fine .* U_Pm_be_fine;
U_I_be_fine  = fresnel_prop_pupil_to_image(U_P_be_fine, D, delta_f, wvl, z_I);
U_S_be_fine  = sqrt(m_S_be)*U_I_be_fine/sqrt(mean(abs(U_I_be_fine(:)).^2));

% Coarse-grid versions for Fig 3 (e,f) phase alignment reference only
U_S_be = U_S_be_fine(1:K:end, 1:K:end);
O      = O_fine;
R_O    = R_O_fine;

[dH_be_ideal, dH_be] = record_holograms_ippsrg(U_S_be_fine, U_R_cell, B, ...
    sigma_r, ell, g_AD, nBits, K);
holoMaxDN_be = max(dH_be(:));
dH_be_ideal = register_pixel_centres(dH_be_ideal, K);
dH_be       = register_pixel_centres(dH_be, K);

dH0_be_ideal = dH_be_ideal - round(mean(dH_be_ideal(:,:,1),'all'));
dH0_be       = dH_be       - round(mean(dH_be(:,:,1),'all'));

USe_be0 = reconstruct_estimate_psrg(dH0_be_ideal/g_AD, ones(N));
USe_be  = reconstruct_estimate_psrg(dH0_be/g_AD,       ones(N));

% Pupil estimate (same path as the point source); truth = bulls-eye pupil field
U_P_be_coarse = U_P_be_fine(1:K:end, 1:K:end);
USe_be0_pad = zeros(Npad); USe_be0_pad(ip,ip) = USe_be0;
USe_be_pad  = zeros(Npad); USe_be_pad(ip,ip)  = USe_be;
tmp0 = ifft2c(USe_be0_pad); U_P_est_be0 = pupilMask .* tmp0(ip,ip);
tmp  = ifft2c(USe_be_pad);  U_P_est_be  = pupilMask .* tmp(ip,ip);
SF_be0 = field_strehl(U_P_be_coarse, U_P_est_be0, pupilMask);
SF_be  = field_strehl(U_P_be_coarse, U_P_est_be,  pupilMask);

eta_be = pixel_eta(U_S_be_fine.*conj(U_R_base), K);
SNR_be = snr_theory(m_S_be, eta_be);

nz_be = ippsrg_noise_stats(USe_be, U_ref_only);
psrg_report('Bulls-eye', U_S_be_fine, U_R_cell, B, K, m_S_be, m_R, holoMaxDN_be, ...
    nBits, eta_be, nz_be, SNR_be, sigma_PS2, kappa, SF_be0, SF_be);

% Fig 3 (a-d): bulls-eye delta = 0 and pi holograms
fig3 = make_figure(figW,figH3);
t3   = tiledlayout(fig3,3,2,'TileSpacing','compact','Padding','compact');
cH_be = max(abs(dH0_be_ideal(:,:,[1 3])),[],'all');   % delta = 0 and pi
ax_ = show_holo(t3,uPix,uPix,dH0_be_ideal(:,:,1),cH_be,labNx,labNy,cbHolo0, '(a)');
add_zoom_inset(ax_,uPix,uPix,dH0_be_ideal(:,:,1),[-48 -32],[-48 -32],[50 124],[50 124]);
ax_ = show_holo(t3,uPix,uPix,dH0_be(:,:,1),      cH_be,labNx,labNy,cbHoloN0,'(b)');
add_zoom_inset(ax_,uPix,uPix,dH0_be(:,:,1),[-48 -32],[-48 -32],[50 124],[50 124]);
USe_be_disp = align_global_phase(quad_phase_remove(U_S_be), USe_be, true(N));
ax_ = show_holo(t3,uPix,uPix,dH0_be_ideal(:,:,3),cH_be,labNx,labNy,cbHoloPi, '(c)');
add_zoom_inset(ax_,uPix,uPix,dH0_be_ideal(:,:,3),[-48 -32],[-48 -32],[50 124],[50 124]);
ax_ = show_holo(t3,uPix,uPix,dH0_be(:,:,3),      cH_be,labNx,labNy,cbHoloNPi,'(d)');
add_zoom_inset(ax_,uPix,uPix,dH0_be(:,:,3),[-48 -32],[-48 -32],[50 124],[50 124]);

% Fig 3 (e,f): truth and estimated image-plane irradiance
Itru = abs(U_I_be_fine).^2;        Itru = Itru/max(Itru(:));
Iest = abs(USe_be_disp).^2;        Iest = Iest/max(Iest(:));
objMask = double(O_fine ~= 0);
labWx='$x/W$'; labWy='$y/W$';
uObj_fine = x2f_vec/D;   % fine grid normalized by W=D
uObj      = x2_vec/D;    % coarse grid normalized by W=D

show_amp  (t3,uObj_fine,uObj_fine,Itru,[],labWx,labWy,cbIrrUS, '(e)');
show_amp  (t3,uObj,     uObj,     Iest,[],labWx,labWy,cbIrrUSe,'(f)');

names = {'IPPSRG_pt','IPPSRG_bulls'};
figs  = [fig1,fig3];
for n = 1:2
    exportgraphics(figs(n),['../outputs/',names{n},'.pdf'], ...
        'ContentType','vector','BackgroundColor','white');
    exportgraphics(figs(n),['../outputs/',names{n},'.png'],'Resolution',600);
end

% ============================================================
% Local functions
% ============================================================

function [dH_ideal, dH] = record_holograms_ippsrg(U_S_fine, U_R_cell, B, ...
    sigma_r, ell, g_AD, nBits, K)
% Record four phase-shifted holograms with pixel binning by K.
    Nfine4 = size(U_S_fine,1);
    Nd     = Nfine4/K;
    dH_ideal = zeros(Nd,Nd,4);
    dH       = zeros(Nd,Nd,4);
    usePoisson = license('test','statistics_toolbox');
    for kk = 1:4
        mH_fine = abs(U_S_fine/sqrt(B) + U_R_cell{kk}).^2;
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

function U_N = reference_only_ippsrg(U_R_cell, sigma_r, ell, g_AD, nBits, K)
    Nfine4   = size(U_R_cell{1},1);
    U_S_zero = zeros(Nfine4);
    [~, dH_n] = record_holograms_ippsrg(U_S_zero, U_R_cell, 1, ...
        sigma_r, ell, g_AD, nBits, K);
    dH_n  = register_pixel_centres(dH_n, K);
    Nd    = size(dH_n,1);
    dH0_n = dH_n - round(mean(dH_n(:,:,1),'all'));
    U_N   = reconstruct_estimate_psrg(dH0_n/g_AD, ones(Nd));
end

function nz = ippsrg_noise_stats(USe_sn, U_N)
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

function psrg_report(name,U_S,U_R_cell,B,K,mS,mR,holoMaxDN,nBits,eta,nz,snrClosed,sigma_PS2,kappa,SF0,SF)
% Standard console block (PRINT_OBJECT_REPORT) for the phase-shifting geometries.
% Noise and signal are per pixel of the signal estimate [pe^2].
    r.mS = mS;  r.mS_meas = mean(pixel_bin(abs(U_S).^2,K),'all');
    r.mR = mR;  r.mR_meas = mean(pixel_bin(abs(U_R_cell{1}).^2,K),'all');
    mx = 0;
    for kk = 1:numel(U_R_cell)
        mH = pixel_bin(abs(U_S/sqrt(B) + U_R_cell{kk}).^2,K);
        mx = max(mx, max(mH(:)));
    end
    r.holoMaxPe = mx;  r.holoMaxDN = holoMaxDN;  r.maxDN = 2^nBits-1;
    r.eta = eta;
    r.noisePred = sigma_PS2;         r.noiseMeas = nz.varNoise;
    r.nEff      = numel(U_S)/K^2;    % noise stats over the full N x N grid
    r.sigPred   = eta*mR*mS/B;       r.sigMeas   = nz.signalPower;
    r.snrClosed = snrClosed;  r.snrNum = nz.SNR;
    r.kappa = kappa;
    r.SF0 = SF0;  r.SF = SF;
    print_object_report(name,r);
end
