function l = gp_lengthscales(th_grid, gp)
% Returns lengthscales 'l' of the GP model.

d = th_grid.dim;
if ~isempty(gp)
    [w,s] = gp_pak(gp);
    l = exp(w(2:(d+1))); % this should work for both GP models
    l = l(:);
else
    l = Inf(d,1);
end
end
