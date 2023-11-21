function [] = suptitle(titl)
% This simple function circumvents an error caused by missing suptitle due
% to Bioinformatics toolbox not being installed.
% Draws no title if suptitle cannot be called for whatever reason.

try
    builtin('suptitle',titl);
catch
    % Do nothing i.e. general title is not printed at all.
    %titl
end
end
