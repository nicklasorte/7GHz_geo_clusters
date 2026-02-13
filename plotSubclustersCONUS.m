function plotSubclustersCONUS(out)
% plotSubclustersCONUS
% Plot subcluster midpoints, limited to continental USA.
tic;
T = out.T;

% Continental USA bounds
latlim = [24 50];
lonlim = [-125 -66];

% Filter to CONUS only
inCONUS = T.latMid >= latlim(1) & T.latMid <= latlim(2) & ...
          T.lonMid >= lonlim(1) & T.lonMid <= lonlim(2);

Tplot = T(inCONUS,:);

gx = geoaxes;
geobasemap(gx,"streets");
hold(gx,"on")

if isempty(Tplot)
    geolimits(gx, latlim, lonlim)
    title(gx,"No midpoints found in CONUS bounds")
    return
end

geoscatter(gx, Tplot.latMid, Tplot.lonMid, 18, Tplot.subclusterID, "filled")

% Force map window to CONUS
geolimits(gx, latlim, lonlim)

% Colormap + colorbar
colormap(gx, turbo(max(Tplot.subclusterID)));
cb = colorbar;                 % <-- key fix
cb.Label.String = "subclusterID";

title(gx, sprintf("Subclusters (CONUS) - %d clusters shown", max(Tplot.subclusterID)))
toc;
end
