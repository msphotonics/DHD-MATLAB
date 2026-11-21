% ============================================================
% script_IPRG_illustrative.m
%
% Illustrative examples for the off-axis IPRG (Section 4.5).
%
% Produces four figures:
%   Figure 1 (3x2): point-source hologram and Fourier plane
%   Figure 2 (2x2): point-source truth vs. estimate, pupil plane
%   Figure 3 (3x2): bulls-eye hologram and Fourier plane
%   Figure 4 (2x2): bulls-eye truth vs. estimate, image plane
%
% Axis conventions:
%   pupil / Fourier planes:  x_1/D,  y_1/D
%   image plane (FPA):       x_2/p,  y_2/p
%
% All functions are in ../functions/.
% Figures and results are saved to ../outputs/.
% ============================================================

clear; close all; clc;
rng(1);

% Run from this script's folder, so ../functions and ../outputs resolve
% wherever MATLAB's current folder happens to be
cd(fileparts(mfilename('fullpath')));
addpath('../functions');

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------

N   = 256;  K = 2;  q_I = 3;

D   = 0.30;          wvl = 1.0e-6;
phi_rms = 1.0;       phaseCorr = D/10;
sigma_r = 100;       ell = 100e3;      m_R = 0.75*ell;
m_S_point = 1;       m_S_bars  = 4;
nBits = 14;          g_AD = (2^nBits-1)/ell;
sigma_q2 = 1/12;

enforceWellDepth = true;
usePoisson = license('test','statistics_toolbox');

holoZoom = 48;       % half-width [px] of point-source hologram zoom
dBmin    = -80;

figW  = 6.5;  figH3 = 6.5;  figH2 = 5.0;  figH1 = 2.7;
fontSize = 10;

exportFigures = true;
saveResults   = true;

set(groot, ...
    'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultColorbarTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',fontSize, ...
    'defaultTextFontSize',fontSize, ...
    'defaultColorbarFontSize',fontSize, ...
    'defaultAxesLineWidth',0.5);

% ------------------------------------------------------------
% Derived quantities
% ------------------------------------------------------------

M = round(q_I*N);
if abs(q_I*N-M)>1e-9 || mod(M,4)~=0
    error('q_I*N must be an integer divisible by 4.');
end
Nfine = K*N;  Mfine = K*M;

delta1 = D/N;  p = delta1;  z_I = q_I*D*p/wvl;
sigma_n2 = m_R + sigma_r^2 + sigma_q2/g_AD^2;

print_run_header('OFF-AXIS IPRG ILLUSTRATIVE EXAMPLES',N,K,'q_I',q_I,D,wvl,'z_I',z_I,p,nBits,g_AD,sigma_n2);

% ------------------------------------------------------------
% Pupil plane
% ------------------------------------------------------------

% Fine pupil grid (Nfine×Nfine at delta_f)
delta_f = delta1/K;
x1f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[x1fGrid,y1fGrid] = meshgrid(x1f_vec,x1f_vec);
pupilMask_fine = double(hypot(x1fGrid,y1fGrid) <= D/2);
phi_fine = make_isoplanatic_phase(Nfine,delta_f,phaseCorr,phi_rms,pupilMask_fine);
P_fine   = pupilMask_fine.*exp(1j*phi_fine);
% Coarse pupil (N×N): alias-free decimation
x1 = x1f_vec(1:K:end);
[x1Grid,y1Grid] = meshgrid(x1,x1);
pupilMask = pupilMask_fine(1:K:end,1:K:end);
P_coarse  = P_fine(1:K:end,1:K:end);
delta1    = delta_f*K;  p = delta1;

% Fine image-plane reference — Eq. (28): spherical wave + tilt
% Grid spacing = delta_f = p/K, sideband at (+M/4, +M/4)
kwave = 2*pi/wvl;
x2f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[X2,Y2] = meshgrid(x2f_vec,x2f_vec);
x_R = D*q_I/4;   % tilt places sideband at +M/4 in Fourier plane
U_R_fine = (1/1j)*exp(1j*kwave*z_I) ...
         .*exp(+1j*kwave/(2*z_I)*(X2.^2+Y2.^2)) ...
         .*exp(+1j*2*pi*x_R*(X2+Y2)/(wvl*z_I));
U_R_fine = sqrt(m_R)*U_R_fine/sqrt(mean(abs(U_R_fine(:)).^2));

% Sideband at (+M/4,+M/4) under ifft2c
sbRow = M/2+1+M/4;  sbCol = M/2+1+M/4;
rowIndex = (sbRow-N/2):(sbRow+N/2-1);
colIndex = (sbCol-N/2):(sbCol+N/2-1);

% Circular sideband window mask (radius N/2), matching the original paper.
% Used for both signal and reference-only noise estimation.
[cwX,cwY] = meshgrid(linspace(-(N-1)/2,(N-1)/2,N),linspace(-(N-1)/2,(N-1)/2,N));
circWin   = hypot(cwX,cwY) <= (N-1)/2;

% Closed-form SNR (Eq. 67 form); eta from PIXEL_ETA on the cross term
snr_theory    = @(mS,eta) eta*(4*q_I^2/pi)*m_R*mS/sigma_n2;
sigma_n2_full = @(mS) mS + sigma_n2;

% Display axes
uPup  = x1/D;
uPix  = (0:N-1)-N/2;
uFour = ((0:M-1)-M/2)/N;   % x_1/D units (= cycles/pixel * q_I)
uObj  = (0:M-1)-M/2;

labPx = '$x_1/D$';  labPy = '$y_1/D$';
labNx = '$x_2/p$';  labNy = '$y_2/p$';

% Colorbar labels
cbAmpUP   = '$|U_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUP   = '$\angle U_\mathrm{P}\;(\mathrm{rad})$';
cbAmpUPe  = '$|\hat{U}_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUPe  = '$\angle\hat{U}_\mathrm{P}\;(\mathrm{rad})$';
cbHolo    = '$d_\mathrm{H}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloN   = '$d_\mathrm{H}^+-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbFourDB  = '$|\widetilde{\rho}_\mathrm{H}^{+}|\;(\mathrm{dB})$';
cbFourPhz = '$\angle\widetilde{\rho}_\mathrm{H}^{+}\;(\mathrm{rad})$';
cbAmpRO   = '$|R_\mathrm{O}|\;(\mathrm{norm.})$';
cbPhzRO   = '$\angle R_\mathrm{O}\;(\mathrm{rad})$';
cbIrrUI   = '$|U_\mathrm{I}|^2\;(\mathrm{norm.})$';
cbIrrUIe  = '$|\hat{U}_\mathrm{I}|^2\;(\mathrm{norm.})$';

% ============================================================
% PART I: POINT-SOURCE OBJECT
% ============================================================

U_P_point = P_coarse;  % N×N truth (decimated from fine grid)

% Fresnel propagate fine pupil to fine image plane (no tilt on pupil)
% U_I_pt_fine: unnormalized, for isoplanatic convolution later
% U_I_pt:      normalized, for hologram recording
U_I_pt_fine = fresnel_prop_pupil_to_image(P_fine, D, delta_f, wvl, z_I);
% m_S = mean over the FULL N x N FPA, as in the other three scripts
U_I_pt      = sqrt(m_S_point/mean(pixel_bin(abs(U_I_pt_fine).^2,K),'all'))*U_I_pt_fine;

[mH_pt,dH_pt_ideal,dH_pt] = record_holograms(U_I_pt,{U_R_fine},K, ...
    sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits);
mH_pt = squeeze(mH_pt);  dH_pt_ideal = squeeze(dH_pt_ideal);  dH_pt = squeeze(dH_pt);
holoMaxDN_pt = max(dH_pt(:));
% Binned values sit at pixel centres; register them to the nominal grid
dH_pt_ideal = register_pixel_centres(dH_pt_ideal,K);
dH_pt       = register_pixel_centres(dH_pt,K);

dH0_pt_ideal = dH_pt_ideal-round(mean(dH_pt_ideal(:)));
dH0_pt       = dH_pt      -round(mean(dH_pt(:)));

[UPe_pt0,F_pt0] = reconstruct_estimate_iprg(dH0_pt_ideal/g_AD,M,rowIndex,colIndex,pupilMask);
[UPe_pt, F_pt ] = reconstruct_estimate_iprg(dH0_pt/g_AD,      M,rowIndex,colIndex,pupilMask);

check_sideband(U_P_point,UPe_pt0,pupilMask,'Point source');

SF_pt0 = field_strehl(U_P_point,UPe_pt0,pupilMask);
SF_pt  = field_strehl(U_P_point,UPe_pt, pupilMask);
eta_pt = pixel_eta(U_I_pt.*conj(U_R_fine),K);
SNR_pt = snr_theory(m_S_point,eta_pt);

% Reference-only noise estimate (signal = 0), matching the trilogy papers.
[~,F_ref_pt] = reconstruct_estimate_iprg( ...
    reference_only_hologram(U_R_fine,K,sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits), ...
    M,rowIndex,colIndex,pupilMask);
nz_pt  = ref_noise_stats(F_pt, F_ref_pt, rowIndex,colIndex,circWin);
nz_pt0 = ref_noise_stats(F_pt0,F_ref_pt, rowIndex,colIndex,circWin);
offaxis_report('Point source',U_I_pt,U_R_fine,K,m_S_point,m_R,mH_pt,holoMaxDN_pt, ...
    nBits,eta_pt,nz_pt,SNR_pt,sigma_n2_full(m_S_point),q_I,N,SF_pt0,SF_pt);

UPe_pt_a = align_global_phase(U_P_point,UPe_pt,pupilMask);

% ------------------------------------------------------------
% Figure 1 (3x2): point-source hologram and Fourier plane
% ------------------------------------------------------------

fig1 = make_figure(figW,figH2);
t1   = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');

zi = centered_indices(N,2*holoZoom);
cH = max(abs([dH0_pt_ideal(zi,zi); dH0_pt(zi,zi)]),[],'all');

ax_ = show_holo(t1,uPix(zi),uPix(zi),dH0_pt_ideal(zi,zi),cH,labNx,labNy,cbHolo, '(a)');
add_zoom_inset(ax_,uPix(zi),uPix(zi),dH0_pt_ideal(zi,zi),[-8 8],[-8 8],[12 46],[12 46]);
ax_ = show_holo(t1,uPix(zi),uPix(zi),dH0_pt(zi,zi),      cH,labNx,labNy,cbHoloN,'(b)');
add_zoom_inset(ax_,uPix(zi),uPix(zi),dH0_pt(zi,zi),[-8 8],[-8 8],[12 46],[12 46]);
show_fourier_db   (t1,uFour,F_pt,q_I,dBmin,cbFourDB, '(c)');
show_fourier_phase(t1,uFour,F_pt,q_I,       cbFourPhz,'(d)');

% ------------------------------------------------------------
% Figure 2 (2x2): point-source truth vs. estimate, pupil plane
% ------------------------------------------------------------

fig2 = make_figure(figW,figH1);
t2   = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

show_phase(t2,uPup,uPup,U_P_point,     pupilMask,labPx,labPy,cbPhzUP, '(a)');
show_phase(t2,uPup,uPup,UPe_pt_a,      pupilMask,labPx,labPy,cbPhzUPe,'(b)');

% ============================================================
% PART II: EXTENDED BULLS-EYE OBJECT
% ============================================================

% Bulls-eye on fine grid — physical propagation (object → pupil → image)
[xOf,yOf] = meshgrid((-Nfine/2:Nfine/2-1)*delta_f);
O_fine   = make_bullseye_target_phys(xOf,yOf,D);
R_O_fine = O_fine.*exp(1j*(-pi+2*pi*rand(Nfine)));
% Step 1: propagate object to pupil (z_P = z_I for unit magnification)
U_Pm_fine    = fresnel_prop_object_to_pupil(R_O_fine, D, delta_f, wvl, z_I);
% Step 2: apply aberrated pupil mask
U_P_be_fine  = P_fine .* U_Pm_fine;
% Step 3: propagate pupil to image plane
U_I_b_raw    = fresnel_prop_pupil_to_image(U_P_be_fine, D, delta_f, wvl, z_I);
% Coarse versions for display (alias-free decimation)
O   = O_fine(1:K:end,1:K:end);
R_O = R_O_fine(1:K:end,1:K:end);
% Pupil truth for S_F: the bulls-eye pupil field (alias-free decimation)
U_P_be_coarse = U_P_be_fine(1:K:end,1:K:end);
U_I_b = sqrt(m_S_bars/mean(pixel_bin(abs(U_I_b_raw).^2,K),'all'))*U_I_b_raw;
% (U_I_b is Nfine×Nfine — ready for record_holograms with pixel_bin)

[mH_b,dH_b_ideal,dH_b] = record_holograms(U_I_b,{U_R_fine},K, ...
    sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits);
mH_b = squeeze(mH_b);  dH_b_ideal = squeeze(dH_b_ideal);  dH_b = squeeze(dH_b);
holoMaxDN_b = max(dH_b(:));
dH_b_ideal = register_pixel_centres(dH_b_ideal,K);
dH_b       = register_pixel_centres(dH_b,K);

dH0_b_ideal = dH_b_ideal-round(mean(dH_b_ideal(:)));
dH0_b       = dH_b      -round(mean(dH_b(:)));

[UPe_b0,F_b0] = reconstruct_estimate_iprg(dH0_b_ideal/g_AD,M,rowIndex,colIndex,pupilMask);
[UPe_b, F_b ] = reconstruct_estimate_iprg(dH0_b/g_AD,      M,rowIndex,colIndex,pupilMask);

check_sideband(U_P_be_coarse,UPe_b0,pupilMask,'Bulls-eye');

SF_b0 = field_strehl(U_P_be_coarse,UPe_b0,pupilMask);
SF_b  = field_strehl(U_P_be_coarse,UPe_b, pupilMask);

eta_b = pixel_eta(U_I_b.*conj(U_R_fine),K);
SNR_b = snr_theory(m_S_bars,eta_b);

% Reuse the same reference-only Fourier plane (same reference, same noise)
nz_b  = ref_noise_stats(F_b, F_ref_pt, rowIndex,colIndex,circWin);
nz_b0 = ref_noise_stats(F_b0,F_ref_pt, rowIndex,colIndex,circWin);
offaxis_report('Bulls-eye',U_I_b,U_R_fine,K,m_S_bars,m_R,mH_b,holoMaxDN_b, ...
    nBits,eta_b,nz_b,SNR_b,sigma_n2_full(m_S_bars),q_I,N,SF_b0,SF_b);

UPe_b_a = align_global_phase(U_P_be_coarse,UPe_b,pupilMask);

% Truth: fine-grid irradiance (before normalization)
Itru = abs(U_I_b_raw).^2; Itru = Itru/max(Itru(:));
% Estimate: zero-pad to q_I*N, fft2c, crop central N×N (matches PP-PSRG)
Npad = M;  % q_I*N
UPe_pad = zeros(Npad);
ip_pad = round(Npad/2-N/2)+(1:N);
UPe_pad(ip_pad,ip_pad) = UPe_b_a;
U_I_est_full = fft2c(UPe_pad);
ic = round(Npad/2-N/2)+(1:N);
U_I_est = U_I_est_full(ic,ic);
Iest = abs(U_I_est).^2; Iest = Iest/max(Iest(:));
objMask = double(O ~= 0);
labWx='$x/W$'; labWy='$y/W$';
uObj_fine   = ((0:Nfine-1)-Nfine/2)*delta_f/D;
uObj_coarse = x1/D;

% ------------------------------------------------------------
% Figure 3 (3x2): bulls-eye hologram and Fourier plane
% ------------------------------------------------------------

fig3 = make_figure(figW,figH2);
t3   = tiledlayout(fig3,2,2,'TileSpacing','compact','Padding','compact');

cHb = max(abs([dH0_b_ideal;dH0_b]),[],'all');

ax_ = show_holo(t3,uPix,uPix,dH0_b_ideal,cHb,labNx,labNy,cbHolo, '(a)');
add_zoom_inset(ax_,uPix,uPix,dH0_b_ideal,[-48 -32],[-48 -32],[50 124],[50 124]);
ax_ = show_holo(t3,uPix,uPix,dH0_b,      cHb,labNx,labNy,cbHoloN,'(b)');
add_zoom_inset(ax_,uPix,uPix,dH0_b,[-48 -32],[-48 -32],[50 124],[50 124]);
show_fourier_db   (t3,uFour,F_b,q_I,dBmin,cbFourDB, '(c)');
show_fourier_phase(t3,uFour,F_b,q_I,       cbFourPhz,'(d)');

% ------------------------------------------------------------
% Figure 4 (2x2): bulls-eye truth vs. estimate, image plane
% ------------------------------------------------------------

fig4 = make_figure(figW,figH1);
t4   = tiledlayout(fig4,1,2,'TileSpacing','compact','Padding','compact');

show_amp  (t4,uObj_fine,  uObj_fine,  Itru,   [],     labWx,labWy,cbIrrUI, '(a)');
show_amp  (t4,uObj_coarse,uObj_coarse,Iest,   [],     labWx,labWy,cbIrrUIe,'(b)');  % N×N crop from padded fft2c

% ------------------------------------------------------------
% Export
% ------------------------------------------------------------

if exportFigures
    names = {'IPRG_pt_hologram','IPRG_pt_comparison', ...
             'IPRG_bars_hologram','IPRG_bars_comparison'};
    figs  = [fig1,fig2,fig3,fig4];
    for n = 1:4
        exportgraphics(figs(n),['../outputs/',names{n},'.pdf'], ...
            'ContentType','vector','BackgroundColor','white');
        exportgraphics(figs(n),['../outputs/',names{n},'.png'],'Resolution',600);
    end
end

% ------------------------------------------------------------
% Save results for Section 7
% ------------------------------------------------------------

if saveResults
    resultIPRG = struct( ...
        'N',N,'K',K,'q_I',q_I,'D',D,'wvl',wvl,'z_I',z_I,'p',p, ...
        'phi_rms',phi_rms,'sigma_r',sigma_r,'ell',ell,'m_R',m_R, ...
        'nBits',nBits,'g_AD',g_AD,'sigma_q2',sigma_q2, ...
        'sigma_n2',sigma_n2, ...
        'pupilMask',pupilMask,'phi_fine',phi_fine);
    resultIPRG.point = struct('m_S',m_S_point,'truth',U_P_point, ...
        'est',UPe_pt,'est0',UPe_pt0,'SF',SF_pt,'SF0',SF_pt0, ...
        'eta',eta_pt,'SNRtheory',SNR_pt,'hologramPe',mH_pt,'hologramDN',dH_pt, ...
        'hologramDNnoiseless',dH_pt_ideal, ...
        'noiseStats',nz_pt);
    resultIPRG.bullseye = struct('m_S',m_S_bars,'object',R_O, ...
        'truth',U_P_be_coarse,'est',UPe_b,'est0',UPe_b0,'SF',SF_b,'SF0',SF_b0, ...
        'SNRtheory',SNR_b,'eta',eta_b, ...
        'imageTruth',Itru,'imageEst',Iest,'hologramPe',mH_b, ...
        'hologramDN',dH_b,'hologramDNnoiseless',dH_b_ideal, ...
        'noiseStats',nz_b);
    save('../outputs/IPRG_illustrative_results.mat','resultIPRG');
end

% ============================================================
% Local functions
% ============================================================

function O = make_bullseye_target_phys(xO, yO, D)
    O = double(hypot(xO,yO) <= D/2) ...
      - double(hypot(xO,yO) <= 3*D/8) ...
      - double(hypot(xO,yO) <= D/4) ...
      + double(hypot(xO,yO) <= D/8);
end

function check_sideband(truth,est,mask,name)
    SFd = field_strehl(truth,est,mask);
    SFc = field_strehl(truth,conj(rot90(est,2)),mask);
    if SFc > SFd
        error('%s: conjugate sideband (S_F=%.4f vs %.4f). Check reference sign.',name,SFc,SFd);
    end
end

function holoRef = reference_only_hologram(U_R_fine,K,sigma_r,usePoisson,ell,enforce,g,nBits)
% Record a hologram with signal = 0 for reference-only noise estimation.
    Nfine  = size(U_R_fine,1);
    U_zero = zeros(Nfine,'like',U_R_fine);
    mH_ref = pixel_bin(abs(U_zero + U_R_fine).^2,K);
    dH_ref = detect_dn(mH_ref,sigma_r,usePoisson,ell,enforce,g,nBits);
    dH_ref = register_pixel_centres(dH_ref,K);
    holoRef = (dH_ref - round(mean(dH_ref(:)))) / g;
end

function nz = ref_noise_stats(F_sn, F_ref, rowIndex, colIndex, circWin)
% Numerical SNR from reference-only noise estimate, matching trilogy papers.
%   F_sn  - Fourier plane of signal+noise hologram
%   F_ref - Fourier plane of reference-only hologram (no signal)
% Both evaluated within the circular sideband window (circWin).
    in   = circWin;
    s_n  = F_sn( rowIndex,colIndex); s_n  = s_n(in);
    n    = F_ref(rowIndex,colIndex); n    = n(in);
    nz.powerSN      = mean(abs(s_n).^2);
    nz.powerN       = mean(abs(n).^2);
    nz.varNoise     = var(n);
    nz.signalPower  = nz.powerSN - nz.powerN;
    nz.SNR          = nz.signalPower / nz.varNoise;
end

function offaxis_report(name,U_S,U_R,K,mS,mR,mH,holoMaxDN,nBits,eta,nz,snrClosed,sig2,q,N,SF0,SF)
% Standard console block (PRINT_OBJECT_REPORT) for the off-axis geometries.
% Noise and signal predictions are in the units of the windowed sideband.
    r.mS = mS;        r.mS_meas = mean(pixel_bin(abs(U_S).^2,K),'all');
    r.mR = mR;        r.mR_meas = mean(pixel_bin(abs(U_R).^2,K),'all');
    r.holoMaxPe = max(mH(:));  r.holoMaxDN = holoMaxDN;  r.maxDN = 2^nBits-1;
    r.eta = eta;
    r.noisePred = sig2/(q^4*N^2);           r.noiseMeas = nz.varNoise;
    r.nEff      = pi*N^2/(4*q^2);
    r.sigPred   = eta*4*mS*mR/(pi*q^2*N^2); r.sigMeas   = nz.signalPower;
    r.snrClosed = snrClosed;  r.snrNum = nz.SNR;
    r.kappa = 1;              % SNR already referred to the pupil support
    r.SF0 = SF0;  r.SF = SF;
    print_object_report(name,r);
end
