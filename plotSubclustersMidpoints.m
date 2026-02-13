function plotSubclustersMidpoints(out)
% plotSubclustersMidpoints
% Plot midpoints colored by subclusterID on a geoaxes map.
%
% Input:
%   out: struct returned by greedySubclustersByRadius (must contain out.T)

T = out.T;

gx = geoaxes;
geobasemap(gx, "streets");
hold(gx, "on")

% Use subclusterID as color data (categorical coloring via numeric mapping).
geoscatter(gx, T.latMid, T.lonMid, 18, T.subclusterID, "filled")

title(gx, sprintf("Subclusters by midpoint (N=%d, clusters=%d)", ...
    height(T), max(T.subclusterID)))

% Improve color variety
colormap(gx, turbo(max(T.subclusterID)));  % turbo gives more distinct colors than default
cb = colorbar(gx);
cb.Label.String = "subclusterID";

end
