function show_phase(t, xv, yv, U, mask, xlab, ylab, cbLab, panel)
% SHOW_PHASE  Display the wrapped phase of a complex field.
%   SHOW_PHASE(T, XV, YV, U, MASK, XLAB, YLAB, CBLAB, PANEL) adds
%   a tile to the tiled layout T showing angle(U) over [-pi, pi],
%   with pixels outside MASK set to NaN.
%
%   Inputs:
%       t      - Tiled layout handle
%       xv, yv - Axis coordinate vectors
%       U      - Complex field array
%       mask   - Binary support mask, or [] for full array
%       xlab   - x-axis label (LaTeX)
%       ylab   - y-axis label (LaTeX)
%       cbLab  - Colorbar label (LaTeX)
%       panel  - Panel letter string, e.g. '(b)'

ph = angle(U);
if ~isempty(mask)
    ph(mask <= 0) = NaN;
end
ax = base_axes(t, xv, yv, ph, xlab, ylab, panel);
clim(ax, [-pi, pi]);
cb            = colorbar(ax);
cb.Ticks      = [-pi, 0, pi];
cb.TickLabels = {'$-\pi$', '$0$', '$\pi$'};
set_cb_label(cb, cbLab);
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
