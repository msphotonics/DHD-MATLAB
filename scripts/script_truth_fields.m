% script_truth_fields.m
%
% Truth fields shared by the illustrative examples in "Digital-Holographic
% Detection" (Appendix C):
%   (a) amplitude and (b) phase of the bulls-eye object-plane complex
%       reflectance R_O, and
%   (c) amplitude and (d) phase of the truth pupil-plane complex field
%       for the point-source object (the generalized pupil function P), and
%   (e) amplitude and (f) phase of the truth image-plane complex field
%       for the bulls-eye object in (a) and (b), i.e., R_O propagated to the
%       pupil, multiplied by P, and propagated to the image plane.
%
% The parameters match the four illustrative scripts. Because each of those
% scripts calls rng(1) and then draws the isoplanatic phase first, the pupil
% phase shown here is the one used in all four. The random phase of R_O is
% drawn later in each script, so panel (b) shows one representative
% realization.

clear; close all; clc;
rng(1);
% Run from this script's folder, so ../functions and ../outputs resolve
% wherever MATLAB's current folder happens to be
cd(fileparts(mfilename('fullpath')));
addpath('../functions');

% ------------------------------------------------------------
% Parameters (match the illustrative scripts)
% ------------------------------------------------------------
N = 256;  K = 2;  q_I = 3;
D = 0.30;  W = D;  wvl = 1.0e-6;
phi_rms = 1.0;  phaseCorr = D/10;

figW = 6.5;  figH3 = 6.5;  fontSize = 10;
exportFigures = true;

set(groot,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultColorbarTickLabelInterpreter','latex', ...
    'defaultAxesFontSize',fontSize,'defaultTextFontSize',fontSize, ...
    'defaultColorbarFontSize',fontSize,'defaultAxesLineWidth',0.5);

% ------------------------------------------------------------
% Grids
% ------------------------------------------------------------
Nfine   = K*N;
delta   = D/N;
delta_f = delta/K;
p       = delta;
z_I     = q_I*D*p/wvl;                     % z_P = z_I (unit magnification)
xf      = (-Nfine/2:Nfine/2-1)*delta_f;
[Xf,Yf] = meshgrid(xf,xf);

% ------------------------------------------------------------
% Truth pupil (same first random draw as the illustrative scripts)
% ------------------------------------------------------------
pupilMask_fine = double(hypot(Xf,Yf) <= D/2);
phi_fine = make_isoplanatic_phase(Nfine,delta_f,phaseCorr,phi_rms,pupilMask_fine);
P_fine   = pupilMask_fine.*exp(1j*phi_fine);
P        = P_fine(1:K:end,1:K:end);            % alias-free decimation
pupilMask = pupilMask_fine(1:K:end,1:K:end);

% ------------------------------------------------------------
% Bulls-eye object-plane complex reflectance
% ------------------------------------------------------------
O_fine   = make_bullseye_target_phys(Xf,Yf,D);
R_O_fine = O_fine.*exp(1j*(-pi+2*pi*rand(Nfine)));
R_O      = R_O_fine(1:K:end,1:K:end);
objMask  = double(R_O ~= 0);

% ------------------------------------------------------------
% Bulls-eye image-plane complex field (object -> pupil -> image),
% the same propagation as in the illustrative scripts
% ------------------------------------------------------------
U_P_b = P_fine.*fresnel_prop_object_to_pupil(R_O_fine,D,delta_f,wvl,z_I);
U_I_b = fresnel_prop_pupil_to_image(U_P_b,D,delta_f,wvl,z_I);

% ------------------------------------------------------------
% Figure
% ------------------------------------------------------------
x   = (-N/2:N/2-1)*delta;
uO  = x/W;   uP = x/D;
labOx = '$x_0/W$';  labOy = '$y_0/W$';
labPx = '$x_1/D$';  labPy = '$y_1/D$';
cbAmpRO = '$|R_\mathrm{O}|$';
cbPhzRO = '$\angle R_\mathrm{O}\;(\mathrm{rad})$';
cbAmpP  = '$|U_\mathrm{P}|\;(\mathrm{norm.})$';
cbPhzP  = '$\angle U_\mathrm{P}\;(\mathrm{rad})$';
labIx = '$x/W$';  labIy = '$y/W$';
uI    = xf/W;                              % fine image-plane grid
cbAmpI  = '$|U_\mathrm{I}|\;(\mathrm{norm.})$';
cbPhzI  = '$\angle U_\mathrm{I}\;(\mathrm{rad})$';

fig1 = make_figure(figW,figH3);
t1   = tiledlayout(fig1,3,2,'TileSpacing','compact','Padding','compact');
show_amp  (t1,uO,uO,abs(R_O),[],       labOx,labOy,cbAmpRO,'(a)');
show_phase(t1,uO,uO,R_O,     objMask,  labOx,labOy,cbPhzRO,'(b)');
show_amp  (t1,uP,uP,abs(P),  pupilMask,labPx,labPy,cbAmpP, '(c)');
show_phase(t1,uP,uP,P,       pupilMask,labPx,labPy,cbPhzP, '(d)');
show_amp  (t1,uI,uI,abs(U_I_b),[],     labIx,labIy,cbAmpI, '(e)');
show_phase(t1,uI,uI,U_I_b,     [],     labIx,labIy,cbPhzI, '(f)');

if exportFigures
    if ~exist('../outputs','dir'), mkdir('../outputs'); end
    exportgraphics(fig1,'../outputs/truth_fields.pdf','ContentType','vector');
    exportgraphics(fig1,'../outputs/truth_fields.png','Resolution',600);
end
fprintf('\n Saved ../outputs/truth_fields.pdf and .png\n');

% ============================================================
function O = make_bullseye_target_phys(xO,yO,D)
    O = double(hypot(xO,yO) <= D/2) ...
      - double(hypot(xO,yO) <= 3*D/8) ...
      - double(hypot(xO,yO) <= D/4) ...
      + double(hypot(xO,yO) <= D/8);
end
