function show_fourier_db(t, uF, F, q, dBmin, cbLab, panel, useSquare, xlab, ylab)
% SHOW_FOURIER_DB  Display Fourier-plane amplitude in dB with overlays.
%   SHOW_FOURIER_DB(T, UF, F, Q, DBMIN, CBLAB, PANEL) adds a tile to
%   the tiled layout T showing 20*log10(|F|/max(|F|)) over [DBMIN, 0]
%   dB, with three overlaid circles marking:
%       - Solid:     positive spatial-heterodyne sideband
%       - Dashed:    adjacent-quadrant noise estimation window
%       - Dash-dot:  two-dimensional autocorrelation support
%
%   Inputs:
%       t      - Tiled layout handle
%       uF     - Fourier-plane axis vector in units of x_1/D
%                (from ((0:M-1)-M/2)/N)
%       F      - Fourier-plane complex array, M x M
%       q      - Sampling quotient (q_I or q_P); controls axis limits
%                and circle positions
%       dBmin  - Lower colorbar limit [dB], e.g. -80
%       cbLab  - Colorbar label (LaTeX)
%       panel  - Panel letter string, e.g. '(e)'
%       useSquare - (optional) square instead of circular sideband markers
%       xlab, ylab - (optional) axis labels; default '$x_1/D$', '$y_1/D$'

if nargin < 9 || isempty(xlab), xlab = '$x_1/D$'; end
if nargin < 10 || isempty(ylab), ylab = '$y_1/D$'; end
A  = 20 * log10(abs(F) / max(abs(F(:))) + eps);
ax = base_axes(t, uF, uF, A, xlab, ylab, panel);
clim(ax, [dBmin, 0]);
set_cb_label(colorbar(ax), cbLab);
if nargin < 8, useSquare = false; end
add_support_overlay(ax, q, useSquare);
end

function add_support_overlay(ax, q, useSquare)
% Overlay for Fourier-plane panels.
% useSquare = true  -> square fiducials (PPRG: rectangular FPA field stop)
% useSquare = false -> circular fiducials (IPRG: circular pupil aperture)
if nargin < 3, useSquare = false; end
hold(ax, 'on');
th = linspace(0, 2*pi, 361);
if useSquare
    hw = 0.5;
    plot(ax,  q/4+hw*[-1 1 1 -1 -1],  q/4+hw*[-1 -1 1 1 -1], 'w-',  'LineWidth',0.75);
    plot(ax,  q/4+hw*[-1 1 1 -1 -1], -q/4+hw*[-1 -1 1 1 -1], 'w--', 'LineWidth',0.75);
else
    plot(ax,  q/4+0.5*cos(th),  q/4+0.5*sin(th), 'w-',  'LineWidth',0.75);
    plot(ax,  q/4+0.5*cos(th), -q/4+0.5*sin(th), 'w--', 'LineWidth',0.75);
end
plot(ax, cos(th), sin(th), 'w-.', 'LineWidth', 0.75);
hold(ax, 'off');
xlim(ax, [-q/2,  q/2]);
ylim(ax, [-q/2,  q/2]);
end

function ax = base_axes(t, xv, yv, A, xlab, ylab, panel)
ax   = nexttile(t);
hImg = imagesc(ax, xv, yv, A);
set(hImg, 'AlphaData', ~isnan(A));
axis(ax, 'image');
set(ax, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'on', ...
    'Layer', 'top', 'Color', 'w');
colormap(ax, turbo(256));
xlabel(ax, xlab);
ylabel(ax, ylab);
text(ax, 0, 1.02, panel, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

function set_cb_label(cb, str)
cb.Label.Interpreter = 'latex';
cb.Label.String      = str;
end
