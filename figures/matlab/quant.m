function q = quant(v, p)
%QUANT  Quantiles (p in 0..1) without the Statistics Toolbox.
v = sort(v(~isnan(v))); n = numel(v);
if n == 0, q = nan(size(p)); return; end
q = interp1(linspace(0, 1, n), v, p);
end
