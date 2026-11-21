function show_fourier_phase(t, uF, F, q, cbLab, panel, useSquare, xlab, ylab)
% SHOW_FOURIER_PHASE  Display Fourier-plane phase with overlays.
%   SHOW_FOURIER_PHASE(T, UF, F, Q, CBLAB, PANEL) adds a tile to
%   the tiled layout T showing the wrapped phase of F over [-pi, pi],
%   with the same three-circle support overlay as SHOW_FOURIER_DB.
%
%   Inputs:
%       t      - Tiled layout handle
%       uF     - Fourier-plane axis vector in units of x_1/D
%       F      - Fourier-plane complex array, M x M
%       q      - Sampling quotient (q_I or q_P)
%       cbLab  - Colorbar label (LaTeX)
%       panel  - Panel letter string, e.g. '(f)'
%       useSquare - (optional) square instead of circular sideband markers
%       xlab, ylab - (optional) axis labels; default '$x_1/D$', '$y_1/D$'

if nargin < 8 || isempty(xlab), xlab = '$x_1/D$'; end
if nargin < 9 || isempty(ylab), ylab = '$y_1/D$'; end
ax            = base_axes(t, uF, uF, angle(F), xlab, ylab, panel);
clim(ax, [-pi, pi]);
cb            = colorbar(ax);
cb.Ticks      = [-pi, 0, pi];
cb.TickLabels = {'$-\pi$', '$0$', '$\pi$'};
set_cb_label(cb, cbLab);
if nargin < 7, useSquare = false; end
add_support_overlay(ax, q, useSquare);
end

function add_support_overlay(ax, q, useSquare)
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
