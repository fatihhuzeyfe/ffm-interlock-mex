function ylab(ax, s)
%YLAB  Bold black y label.
S = ffm_style();
ylabel(ax, s, 'FontName', S.font, 'FontSize', S.fs, 'FontWeight', 'bold', 'Color', 'k');
end
