function show_holo(t, xv, yv, H, c, xlab, ylab, cbLab, panel)
% SHOW_HOLO  Display a zero-mean hologram with symmetric DN scaling.
%   SHOW_HOLO(T, XV, YV, H, C, XLAB, YLAB, CBLAB, PANEL) adds a tile
%   to the tiled layout T showing the integer-DN hologram H with color
%   limits [-C, C] and colorbar ticks at power-of-2 values.
%
%   Inputs:
%       t      - Tiled layout handle
%       xv, yv - Axis coordinate vectors (pixels)
%       H      - Zero-mean hologram array [integer DN]
%       c      - Symmetric color limit [DN]; typically
%                max(abs(H(:)), [], 'all')
%       xlab   - x-axis label (LaTeX)
%       ylab   - y-axis label (LaTeX)
%       cbLab  - Colorbar label (LaTeX)
%       panel  - Panel letter string, e.g. '(c)'

ax = base_axes(t, xv, yv, H, xlab, ylab, panel);
clim(ax, [-c, c]);
k             = floor(log2(max(c, 2)));
cb            = colorbar(ax);
cb.Ticks      = [-2^k, -2^(k-1), 0, 2^(k-1), 2^k];
cb.TickLabels = {sprintf('$-2^{%d}$', k), sprintf('$-2^{%d}$', k-1), ...
                 '$0$', sprintf('$2^{%d}$', k-1), sprintf('$2^{%d}$', k)};
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
