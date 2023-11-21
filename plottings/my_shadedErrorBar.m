function s = my_shadedErrorBar(x,y,lb,ub,col,lw)

if nargin < 5
    col = 'k';
end
if nargin < 6
    lw = 0.5; % linewidth, matlab default 0.5
end

% convert input vectors so that column vectors are also valid inputs
x = x(:)';
y = y(:)';
ub(~isfinite(ub)) = realmax; % ad-hoc fix for some plotting issues caused by Inf-values
lb(~isfinite(lb)) = -realmax;
errBar = [ub(:)'-y;y-lb(:)'];
s = shadedErrorBar(x,y,errBar,'lineprops',{'-','Color',col},'patchSaturation',0.1);
s.mainLine.LineWidth = lw;
box on;
end
