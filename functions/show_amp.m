function show_amp(t, xv, yv, A, mask, xlab, ylab, cbLab, panel)
% SHOW_AMP  Display a normalized amplitude or irradiance panel.
%   SHOW_AMP(T, XV, YV, A, MASK, XLAB, YLAB, CBLAB, PANEL) adds a
%   tile to the tiled layout T showing the array A normalized to its
%   maximum.  Pixels outside MASK are set to NaN (transparent).
%   Pass MASK = [] to display the full array.
%
%   Inputs:
%       t      - Tiled layout handle (from tiledlayout)
%       xv, yv - Axis coordinate vectors for imagesc
%       A      - Real non-negative array to display
%       mask   - Binary support mask (1 = inside, 0 = outside), or []
%       xlab   - x-axis label string (LaTeX)
%       ylab   - y-axis label string (LaTeX)
%       cbLab  - Colorbar label string (LaTeX)
%       panel  - Panel letter string, e.g. '(a)'

if ~isempty(mask)
    A(mask <= 0) = NaN;
end
A  = A / max(A(:), [], 'omitnan');
ax = base_axes(t, xv, yv, A, xlab, ylab, panel);
clim(ax, [0, 1]);
set_cb_label(colorbar(ax), cbLab);
end

% ------------------------------------------------------------------
% Shared sub-functions (also used by show_phase, show_holo, etc.)
% ------------------------------------------------------------------

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
