function h = txt(ax, x, y, s, varargin)
%TXT  Annotation text in data coordinates: black, normal weight, figure font.
S = ffm_style();
h = text(ax, x, y, s, 'FontName', S.font, 'FontSize', S.fsSmall, 'Color', 'k', varargin{:});
end
