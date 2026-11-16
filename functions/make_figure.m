function fig = make_figure(w, h)
% MAKE_FIGURE  Create a publication-sized figure window.
%   FIG = MAKE_FIGURE(W, H) creates a white-background figure
%   window of width W inches and height H inches, with paper
%   properties set to match so that EXPORTGRAPHICS produces
%   a correctly-sized PDF or PNG without whitespace cropping.
%
%   Inputs:
%       w  - Figure width  [inches]
%       h  - Figure height [inches]
%
%   Output:
%       fig - Handle to the new figure

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position',      [1, 1, w, h], ...
    'PaperUnits',    'inches', ...
    'PaperSize',     [w, h], ...
    'PaperPosition', [0, 0, w, h]);
end
