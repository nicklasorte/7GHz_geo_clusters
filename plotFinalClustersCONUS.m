function plotFinalClustersCONUS(outFinal, filename, dpi)
tic;
% plotFinalClustersCONUS
%
% Inputs
%   outFinal : struct returned by buildFinalClustersGreedy
%   filename : string or char (e.g., "final_clusters.jpg")
%   dpi      : resolution (e.g., 300). Default = 300.
%
% Example:
%   plotFinalClustersCONUS(outFinal, "clusters.jpg", 300)

if nargin < 2 || isempty(filename)
    filename = "final_clusters_CONUS.jpg";
end

if nargin < 3 || isempty(dpi)
    dpi = 300;
end

T = outFinal.T;

% Continental USA bounds
latlim = [24 50];
lonlim = [-125 -66];

% Filter to CONUS
inCONUS = T.latMid >= latlim(1) & T.latMid <= latlim(2) & ...
          T.lonMid >= lonlim(1) & T.lonMid <= lonlim(2);

Tplot = T(inCONUS,:);

fig = figure('Color','w');
gx = geoaxes(fig);
geobasemap(gx,"streets");
hold(gx,"on")

if isempty(Tplot)
    geolimits(gx, latlim, lonlim)
    title(gx,"No points in Continental USA")
    exportgraphics(fig, filename, "Resolution", dpi)
    return
end

% Plot clusters
geoscatter(gx, Tplot.latMid, Tplot.lonMid, ...
           18, Tplot.finalClusterID, "filled")

% Highlight centers
cent = Tplot.isCenter;
geoscatter(gx, Tplot.latMid(cent), Tplot.lonMid(cent), ...
           60, "k", "o", "LineWidth", 1.5)

% Set bounds
geolimits(gx, latlim, lonlim)

nClusters = max(Tplot.finalClusterID);
colormap(gx, turbo(nClusters))
cb = colorbar;
cb.Label.String = "finalClusterID";

title(gx, sprintf("Final Clusters (CONUS) — %d clusters", nClusters))
fig.Position = [100 100 1200 900];
pause(0.1)
% --- Save as high-quality JPG ---
exportgraphics(fig, filename, "Resolution", dpi)

fprintf("Saved map to: %s\n", filename);
toc;
end
