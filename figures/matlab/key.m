function key(ax, entries, varargin)
%KEY  Legend from a cell array of {plot-arguments, label} rows, drawn with invisible dummy lines,
%   black text, no box. Example: key(ax, {{'-o', 'Color', c}, 'PLA/PETG'; ...}, 'Location', 'northwest')
S = ffm_style();
h = [];
for i = 1:size(entries, 1)
    a = entries{i, 1};
    h = [h; plot(ax, NaN, NaN, a{:})];
end
lg = legend(ax, h, entries(:, 2), 'FontName', S.font, 'FontSize', S.fsSmall, 'TextColor', 'k', varargin{:});
set(lg, 'Box', 'off');
end
