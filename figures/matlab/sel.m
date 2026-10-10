function idx = sel(col, value)
%SEL  Logical index of rows whose text column equals value.
idx = strcmp(col, value);
end
