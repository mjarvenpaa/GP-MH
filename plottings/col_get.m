function col = col_get(ind)
% Returns some nice colors - useful for plotting and better than MATLAB
% default colors. 
%
% 'ind' must be an integer between 1 and 9. 
% 93 69 12 yellow
% 0 45 74 blue
% 30 75 93 blue2
% 64 8 18 red
% 85 33 10 red2
% 0 0 0 black 
% 50 50 50 gray
% 49 18 56 purple
% 100 0 100 purple2

cs = [93 69 12;0 45 74;30 75 93;64 8 18;85 33 10;0 0 0;50 50 50;49 18 56;100 0 100]/100;
if ~isnumeric(ind) || numel(ind)~=1 || ~any(ind==1:9)
    error('Incorrect input to get_col.');
end
col = cs(ind,:);
end

