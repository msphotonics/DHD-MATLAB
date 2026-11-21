% ============================================================
% script_PPRG_illustrative.m
%
% Illustrative examples for the off-axis PPRG (Section 5.5).
%
% Produces four figures:
%   Figure 1 (3x2): point-source hologram and Fourier plane
%   Figure 2 (2x2): point-source truth vs. estimate, pupil plane
%   Figure 3 (3x2): bulls-eye hologram and Fourier plane
%   Figure 4 (2x2): bulls-eye reflectance, truth image, estimated image
%
% Axis conventions:
%   pupil plane (FPA):   x_1/p,  y_1/p  (Figs 1,3 a-b)
%   Fourier plane:       x_1/D,  y_1/D  (Figs 1,3 c-f)
%   pupil plane:         x_1/D,  y_1/D  (Fig 2)
%   image plane:         x/W,    y/W    (Fig 4)
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

N = 256;  K = 2;  q_P = 3;

D = 0.30;          wvl = 1.0e-6;
phi_rms = 1.0;     phaseCorr = D/10;
sigma_r = 100;     ell = 100e3;      m_R = 0.75*ell;
m_S_point = 1;     m_S_bars = 4;
nBits = 14;        g_AD = (2^nBits-1)/ell;  sigma_q2 = 1/12;

enforceWellDepth = true;
usePoisson = license('test','statistics_toolbox');

dBmin = -80;
figW = 6.5;  figH3 = 6.5;  figH2 = 5.0;  figH1 = 2.7;  fontSize = 10;
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

M = round(q_P*N);
if abs(q_P*N-M)>1e-9 || mod(M,4)~=0
    error('q_P*N must be an integer divisible by 4.');
end
Nfine = K*N;

delta1 = D/N;  p = delta1;  W = D;
z_P = q_P*W*p/wvl;
sigma_n2 = m_R + sigma_r^2 + sigma_q2/g_AD^2;

print_run_header('OFF-AXIS PPRG ILLUSTRATIVE EXAMPLES',N,K,'q_P',q_P,D,wvl,'z_P',z_P,p,nBits,g_AD,sigma_n2);

% ------------------------------------------------------------
% Fine pupil-plane grid (Nfine x Nfine at delta_f)
% ------------------------------------------------------------

delta_f = delta1/K;
x1f_vec = (-Nfine/2:Nfine/2-1)*delta_f;
[x1fGrid,y1fGrid] = meshgrid(x1f_vec,x1f_vec);
pupilMask_fine = double(hypot(x1fGrid,y1fGrid) <= D/2);
phi_fine = make_isoplanatic_phase(Nfine,delta_f,phaseCorr,phi_rms,pupilMask_fine);
P_fine   = pupilMask_fine.*exp(1j*phi_fine);

% Coarse pupil (N x N): alias-free decimation
x1 = x1f_vec(1:K:end);
[x1Grid,y1Grid] = meshgrid(x1,x1);
pupilMask = pupilMask_fine(1:K:end,1:K:end);
P_coarse  = P_fine(1:K:end,1:K:end);

% ------------------------------------------------------------
% Reference wave on fine grid: Eq. (49) tilt-only plane wave
% x_R = D*q_P/4 places sideband at (+M/4,+M/4)
% ------------------------------------------------------------

x_R = D*q_P/4;
[X1f,Y1f] = meshgrid(x1f_vec,x1f_vec);
U_R_fine = exp(-1j*2*pi*x_R*(X1f+Y1f)/(wvl*z_P));
U_R_fine = sqrt(m_R)*U_R_fine/sqrt(mean(abs(U_R_fine(:)).^2));

% Sideband location: N×N block centered at (+M/4,+M/4) in M×M plane
sbRow = M/2+1+M/4;  sbCol = M/2+1+M/4;
rowIndex = (sbRow-N/2):(sbRow+N/2-1);
colIndex = (sbCol-N/2):(sbCol+N/2-1);

% Circular sideband window (radius N/2) for noise stats
[cwX,cwY] = meshgrid(linspace(-(N-1)/2,(N-1)/2,N),linspace(-(N-1)/2,(N-1)/2,N));
circWin = hypot(cwX,cwY) <= (N-1)/2;

% Closed-form SNR (Eq. 67); eta from PIXEL_ETA on the cross term
snr_theory    = @(mS,eta) eta*(4*q_P^2/pi)*m_R*mS/sigma_n2;
sigma_n2_full = @(mS) mS + sigma_n2;

% Display axes
uPup  = x1/D;
uPix  = (0:N-1)-N/2;
uFour = ((0:M-1)-M/2)/N;        % x_1/D for Fourier plane
uIe   = ((0:N-1)-N/2)/N;       % x_2/W for the N x N image-plane estimate (spacing W/N)

labPx='$x_1/D$'; labPy='$y_1/D$';
labNx='$x_1/p$'; labNy='$y_1/p$';
labIx='$x_2/W$'; labIy='$y_2/W$';
labWx='$x/W$';   labWy='$y/W$';

% Hologram zoom insets (pixel units): a 16 x 16-pixel region inside the
% pupil is magnified in the upper-right corner to show the off-axis fringes
zoomX  = [-48 -32];  zoomY  = [-48 -32];
insetX = [ 50 124];  insetY = [ 50 124];

cbAmpUP  = '$|U_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUP  = '$\angle U_\mathrm{P}\;(\mathrm{rad})$';
cbAmpUPe = '$|\hat{U}_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzUPe = '$\angle\hat{U}_\mathrm{P}\;(\mathrm{rad})$';
cbHolo   = '$d_\mathrm{H}-\langle d_\mathrm{H}\rangle\;(\mathrm{DN})$';
cbHoloN  = '$d_\mathrm{H}^+-\langle d_\mathrm{H}^+\rangle\;(\mathrm{DN})$';
cbFourDB = '$|\widetilde{\rho}_\mathrm{H}^{+}|\;(\mathrm{dB})$';
cbFourPhz= '$\angle\widetilde{\rho}_\mathrm{H}^{+}\;(\mathrm{rad})$';
cbAmpUIe = '$|\hat{U}_\mathrm{I}|\;(\mathrm{norm.})$';
cbPhzUIe = '$\angle\hat{U}_\mathrm{I}\;(\mathrm{rad})$';
cbAmpRO  = '$|R_\mathrm{O}|\;(\mathrm{norm.})$';
cbPhzRO  = '$\angle R_\mathrm{O}\;(\mathrm{rad})$';
cbIrrUI  = '$|U_\mathrm{I}|^2\;(\mathrm{norm.})$';
cbIrrUIe = '$|\hat{U}_\mathrm{I}|^2\;(\mathrm{norm.})$';

% Reference-only noise estimate (N×N sideband, no signal)
[UIe_ref,~] = reconstruct_estimate_pprg( ...
    reference_only_hologram(U_R_fine,K,sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits), ...
    M,rowIndex,colIndex,pupilMask);

% ============================================================
% PART I: POINT-SOURCE OBJECT
% ============================================================

% Fine pupil field (FPA is in pupil plane — no propagation needed)
% m_S = mean over the FULL N x N FPA (incl. dark pixels outside pupil),
% as assumed by alpha_P = 4 q_P^2/pi and as in the 2018 script.
U_S_pt_fine = sqrt(m_S_point/mean(pixel_bin(abs(P_fine).^2,K),'all'))*P_fine;

[mH_pt,dH_pt_ideal,dH_pt] = record_holograms(U_S_pt_fine,{U_R_fine},K, ...
    sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits);
mH_pt = squeeze(mH_pt);  dH_pt_ideal = squeeze(dH_pt_ideal);  dH_pt = squeeze(dH_pt);
holoMaxDN_pt = max(dH_pt(:));
% Binned values sit at pixel centres; register them to the nominal grid
dH_pt_ideal = register_pixel_centres(dH_pt_ideal,K);
dH_pt       = register_pixel_centres(dH_pt,K);

dH0_pt_ideal = dH_pt_ideal - round(mean(dH_pt_ideal(:)));
dH0_pt       = dH_pt       - round(mean(dH_pt(:)));

% Reconstruct: N×N image-plane estimate from sideband
[UIe_pt0,F_pt0] = reconstruct_estimate_pprg(dH0_pt_ideal/g_AD,M,rowIndex,colIndex,pupilMask);
[UIe_pt, F_pt ] = reconstruct_estimate_pprg(dH0_pt/g_AD,      M,rowIndex,colIndex,pupilMask);

% Pupil estimate: zero-pad UIe (N×N) to q_P*N, ifft2c, crop central N×N
% Same fix as IP-PSRG: UIe at spacing q_P*delta -> ifft2c -> spacing delta
Npad = q_P*N;
ip_e = centered_indices(Npad,N);   % embed N×N into Npad×Npad
ip_n = centered_indices(Npad,N);   % crop N×N from Npad×Npad (same here)
UIe_pt0_pad = zeros(Npad);  UIe_pt0_pad(ip_e,ip_e) = UIe_pt0;
UIe_pt_pad  = zeros(Npad);  UIe_pt_pad(ip_e,ip_e)  = UIe_pt;
tmp0 = ifft2c(UIe_pt0_pad);  U_P_est_pt0 = pupilMask .* tmp0(ip_n,ip_n);
tmp  = ifft2c(UIe_pt_pad);   U_P_est_pt  = pupilMask .* tmp(ip_n,ip_n);

SF_pt0 = field_strehl(P_coarse,U_P_est_pt0,pupilMask);
SF_pt  = field_strehl(P_coarse,U_P_est_pt, pupilMask);
eta_pt = pixel_eta(U_S_pt_fine.*conj(U_R_fine),K);
SNR_pt = snr_theory(m_S_point,eta_pt);

nz_pt = ref_noise_stats(UIe_pt,UIe_ref,circWin,M);
offaxis_report('Point source',U_S_pt_fine,U_R_fine,K,m_S_point,m_R,mH_pt,holoMaxDN_pt, ...
    nBits,eta_pt,nz_pt,SNR_pt,sigma_n2_full(m_S_point),q_P,N,SF_pt0,SF_pt);

U_P_est_pt_a = align_global_phase(P_coarse,U_P_est_pt,pupilMask);

% Figure 1 (3x2): hologram and Fourier plane
fig1 = make_figure(figW,figH2);
t1   = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');
cH   = max(abs([dH0_pt_ideal;dH0_pt]),[],'all');
ax1a = show_holo(t1,uPix,uPix,dH0_pt_ideal,cH,labNx,labNy,cbHolo, '(a)');
add_zoom_inset(ax1a,uPix,uPix,dH0_pt_ideal,zoomX,zoomY,insetX,insetY);
ax1b = show_holo(t1,uPix,uPix,dH0_pt,      cH,labNx,labNy,cbHoloN,'(b)');
add_zoom_inset(ax1b,uPix,uPix,dH0_pt,zoomX,zoomY,insetX,insetY);
show_fourier_db   (t1,uFour,F_pt,q_P,dBmin,cbFourDB, '(c)',true,labIx,labIy);
show_fourier_phase(t1,uFour,F_pt,q_P,       cbFourPhz,'(d)',true,labIx,labIy);

% Figure 2 (2x2): pupil truth vs. estimate
fig2 = make_figure(figW,figH1);
t2   = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');
show_phase(t2,uPup,uPup,P_coarse,           pupilMask,labPx,labPy,cbPhzUP, '(a)');
show_phase(t2,uPup,uPup,U_P_est_pt_a,       pupilMask,labPx,labPy,cbPhzUPe,'(b)');

% ============================================================
% PART II: EXTENDED BULLS-EYE OBJECT
% ============================================================

% Object on fine grid
[xOf,yOf] = meshgrid((-Nfine/2:Nfine/2-1)*delta_f);
O_fine   = make_bullseye_target_phys(xOf,yOf,D);
R_O_fine = O_fine.*exp(1j*(-pi+2*pi*rand(Nfine)));

% Physical propagation: object -> pupil (collimate) -> pixel bin
U_Pm_fine   = fresnel_prop_object_to_pupil(R_O_fine,D,delta_f,wvl,z_P);
U_P_be_fine = P_fine .* U_Pm_fine;
U_S_be_fine = sqrt(m_S_bars/mean(pixel_bin(abs(U_P_be_fine).^2,K),'all'))*U_P_be_fine;

% Coarse versions for display
O   = O_fine(1:K:end,1:K:end);
R_O = R_O_fine(1:K:end,1:K:end);
U_P_be_coarse = U_P_be_fine(1:K:end,1:K:end);   % pupil truth for S_F

[mH_b,dH_b_ideal,dH_b] = record_holograms(U_S_be_fine,{U_R_fine},K, ...
    sigma_r,usePoisson,ell,enforceWellDepth,g_AD,nBits);
mH_b = squeeze(mH_b);  dH_b_ideal = squeeze(dH_b_ideal);  dH_b = squeeze(dH_b);
holoMaxDN_b = max(dH_b(:));
dH_b_ideal = register_pixel_centres(dH_b_ideal,K);
dH_b       = register_pixel_centres(dH_b,K);

dH0_b_ideal = dH_b_ideal - round(mean(dH_b_ideal(:)));
dH0_b       = dH_b       - round(mean(dH_b(:)));

[UIe_b0,F_b0] = reconstruct_estimate_pprg(dH0_b_ideal/g_AD,M,rowIndex,colIndex,pupilMask);
[UIe_b, F_b ] = reconstruct_estimate_pprg(dH0_b/g_AD,      M,rowIndex,colIndex,pupilMask);

% Pupil estimate (same path as the point source)
UIe_b0_pad = zeros(Npad);  UIe_b0_pad(ip_e,ip_e) = UIe_b0;
UIe_b_pad  = zeros(Npad);  UIe_b_pad(ip_e,ip_e)  = UIe_b;
tmp0 = ifft2c(UIe_b0_pad);  U_P_est_b0 = pupilMask .* tmp0(ip_n,ip_n);
tmp  = ifft2c(UIe_b_pad);   U_P_est_b  = pupilMask .* tmp(ip_n,ip_n);
SF_b0 = field_strehl(U_P_be_coarse,U_P_est_b0,pupilMask);
SF_b  = field_strehl(U_P_be_coarse,U_P_est_b, pupilMask);

eta_b = pixel_eta(U_S_be_fine.*conj(U_R_fine),K);
SNR_b = snr_theory(m_S_bars,eta_b);

nz_b = ref_noise_stats(UIe_b,UIe_ref,circWin,M);
offaxis_report('Bulls-eye',U_S_be_fine,U_R_fine,K,m_S_bars,m_R,mH_b,holoMaxDN_b, ...
    nBits,eta_b,nz_b,SNR_b,sigma_n2_full(m_S_bars),q_P,N,SF_b0,SF_b);

% Fig 4(c): truth image — Fresnel prop fine pupil to fine image plane
U_I_tru_fine = fresnel_prop_pupil_to_image(U_P_be_fine,D,delta_f,wvl,z_P);
Itru = abs(U_I_tru_fine).^2;  Itru = Itru/max(Itru(:));
uObj_fine = (-Nfine/2:Nfine/2-1)*delta_f/D;

% Fig 4(d): estimated image — |UIe_b|^2 (UIe_b is N×N image-plane estimate)
Iest = abs(UIe_b).^2;  Iest = Iest/max(Iest(:));
uObj_est = uIe;   % same N×N axis as Fig 1/3 (e,f)

objMask = double(O ~= 0);
uObj_coarse = x1/D;

% Figure 3 (3x2): bulls-eye hologram and Fourier plane
fig3 = make_figure(figW,figH2);
t3   = tiledlayout(fig3,2,2,'TileSpacing','compact','Padding','compact');
cHb  = max(abs([dH0_b_ideal;dH0_b]),[],'all');
ax3a = show_holo(t3,uPix,uPix,dH0_b_ideal,cHb,labNx,labNy,cbHolo, '(a)');
add_zoom_inset(ax3a,uPix,uPix,dH0_b_ideal,zoomX,zoomY,insetX,insetY);
ax3b = show_holo(t3,uPix,uPix,dH0_b,      cHb,labNx,labNy,cbHoloN,'(b)');
add_zoom_inset(ax3b,uPix,uPix,dH0_b,zoomX,zoomY,insetX,insetY);
show_fourier_db   (t3,uFour,F_b,q_P,dBmin,cbFourDB, '(c)',true,labIx,labIy);
show_fourier_phase(t3,uFour,F_b,q_P,       cbFourPhz,'(d)',true,labIx,labIy);

% Figure 4 (2x2): object, truth image, estimated image
fig4 = make_figure(figW,figH1);
t4   = tiledlayout(fig4,1,2,'TileSpacing','compact','Padding','compact');
show_amp  (t4,uObj_fine,  uObj_fine,  Itru,    [],      labWx,labWy,cbIrrUI, '(a)');
show_amp  (t4,uObj_est,   uObj_est,   Iest,    [],      labWx,labWy,cbIrrUIe,'(b)');

% ------------------------------------------------------------
% Export
% ------------------------------------------------------------

if exportFigures
    names = {'PPRG_pt_hologram','PPRG_pt_comparison', ...
             'PPRG_bulls_hologram','PPRG_bulls_comparison'};
    figs  = [fig1,fig2,fig3,fig4];
    for n = 1:4
        exportgraphics(figs(n),['../outputs/',names{n},'.pdf'], ...
            'ContentType','vector','BackgroundColor','white');
        exportgraphics(figs(n),['../outputs/',names{n},'.png'],'Resolution',600);
    end
end

% ------------------------------------------------------------
% Save results
% ------------------------------------------------------------

if saveResults
    resultPPRG = struct( ...
        'N',N,'K',K,'q_P',q_P,'D',D,'W',W,'wvl',wvl,'z_P',z_P,'p',p, ...
        'phi_rms',phi_rms,'sigma_r',sigma_r,'ell',ell,'m_R',m_R, ...
        'nBits',nBits,'g_AD',g_AD,'sigma_q2',sigma_q2, ...
        'sigma_n2',sigma_n2, ...
        'pupilMask',pupilMask,'phi_fine',phi_fine);
    resultPPRG.point = struct('m_S',m_S_point,'truth',P_coarse, ...
        'est',U_P_est_pt,'est0',U_P_est_pt0,'SF',SF_pt,'SF0',SF_pt0, ...
        'eta',eta_pt,'SNRtheory',SNR_pt,'hologramPe',mH_pt,'hologramDN',dH_pt, ...
        'hologramDNnoiseless',dH_pt_ideal,'noiseStats',nz_pt);
    resultPPRG.bullseye = struct('m_S',m_S_bars,'object',R_O, ...
        'truth',U_P_be_coarse,'est',U_P_est_b,'est0',U_P_est_b0, ...
        'SF',SF_b,'SF0',SF_b0,'SNRtheory',SNR_b,'eta',eta_b, ...
        'imageTruth',Itru,'imageEst',Iest,'hologramPe',mH_b, ...
        'hologramDN',dH_b,'hologramDNnoiseless',dH_b_ideal,'noiseStats',nz_b);
    save('../outputs/PPRG_illustrative_results.mat','resultPPRG');
end

% ============================================================
% Local functions
% ============================================================

function O = make_bullseye_target_phys(xO,yO,D)
    O = double(hypot(xO,yO) <= D/2) ...
      - double(hypot(xO,yO) <= 3*D/8) ...
      - double(hypot(xO,yO) <= D/4) ...
      + double(hypot(xO,yO) <= D/8);
end

function holoRef = reference_only_hologram(U_R_fine,K,sigma_r,usePoisson,ell,enforce,g,nBits)
    mH_ref = pixel_bin(abs(U_R_fine).^2,K);
    dH_ref = detect_dn(mH_ref,sigma_r,usePoisson,ell,enforce,g,nBits);
    dH_ref = register_pixel_centres(dH_ref,K);
    holoRef = (dH_ref - round(mean(dH_ref(:)))) / g;
end

function nz = ref_noise_stats(UIe_sn,UIe_ref,circWin,Mnorm)
    % Average over the image support (circular window of diameter W -> N px),
    % matching the image-average operator <.>_I in Eq. (66)-(67).
    s_n = UIe_sn(circWin)/Mnorm^2;
    n   = UIe_ref(circWin)/Mnorm^2;
    nz.powerSN     = mean(abs(s_n).^2);
    nz.powerN      = mean(abs(n).^2);
    nz.varNoise    = var(n);
    nz.signalPower = nz.powerSN - nz.powerN;
    nz.SNR         = nz.signalPower / nz.varNoise;
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
