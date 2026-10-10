function save_figure(fig, name)
%SAVE_FIGURE  EPS (vector, editable text) and PDF for the journal, 600 dpi PNG for checking.
S = ffm_style();
base = fullfile(S.out, name);
set(fig, 'Renderer', 'painters');
print(fig, [base '.eps'], '-depsc', '-painters');
print(fig, [base '.pdf'], '-dpdf', '-painters');
print(fig, [base '.png'], '-dpng', '-r600');
fprintf('saved %s (.eps .pdf .png)\n', base);
end
