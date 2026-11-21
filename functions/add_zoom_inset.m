function h = add_zoom_inset(ax, xv, yv, A, xZoom, yZoom, xInset, yInset)
% ADD_ZOOM_INSET  Magnified inset of an image region, drawn in the same axes.
%   H = ADD_ZOOM_INSET(AX, XV, YV, A, XZOOM, YZOOM, XINSET, YINSET) marks
%   the region XZOOM x YZOOM (data units) of the image A, displayed in AX
%   with axis vectors XV and YV, and redraws that region magnified over the
%   rectangle XINSET x YINSET (data units) of the same axes. Because the
%   inset is an image object in AX, it shares the colormap and color limits
%   and needs no separate axes, so it exports exactly as displayed.
%
%   Inputs:
%       ax      - Axes that already show A (e.g., from SHOW_HOLO)
%       xv, yv  - Axis vectors used to display A
%       A       - Displayed image (same array passed to imagesc)
%       xZoom   - [xmin xmax] of the region to magnify
%       yZoom   - [ymin ymax] of the region to magnify
%       xInset  - [xmin xmax] where the magnified region is drawn
%       yInset  - [ymin ymax] where the magnified region is drawn
%
%   Output:
%       h       - Struct of graphics handles (box, inset image, frame, lines)

ix = find(xv >= xZoom(1) & xv <= xZoom(2));
iy = find(yv >= yZoom(1) & yv <= yZoom(2));
sub = A(iy, ix);

hold(ax, 'on');
% Box around the magnified region
h.box = rectangle(ax, 'Position', [xZoom(1), yZoom(1), diff(xZoom), diff(yZoom)], ...
    'EdgeColor', 'w', 'LineWidth', 1);
% Magnified copy (pixel centres spread evenly over the inset rectangle)
dx = diff(xInset)/numel(ix);  dy = diff(yInset)/numel(iy);
h.img = image(ax, 'XData', [xInset(1)+dx/2, xInset(2)-dx/2], ...
                  'YData', [yInset(1)+dy/2, yInset(2)-dy/2], ...
                  'CData', sub, 'CDataMapping', 'scaled');
% Frame around the inset
h.frame = rectangle(ax, 'Position', [xInset(1), yInset(1), diff(xInset), diff(yInset)], ...
    'EdgeColor', 'w', 'LineWidth', 1.25);
% Connector lines between corresponding corners
h.lines = [plot(ax, [xZoom(1) xInset(1)], [yZoom(2) yInset(2)], 'w-', 'LineWidth', 0.75), ...
           plot(ax, [xZoom(2) xInset(2)], [yZoom(1) yInset(1)], 'w-', 'LineWidth', 0.75)];
hold(ax, 'off');
end
