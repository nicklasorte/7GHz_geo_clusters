% -------------------------------------------------------------------------
function [latMid, lonMid] = geodesicMidpoints(lat1, lon1, lat2, lon2)
% Great-circle midpoint using Mapping Toolbox.
[arcDeg, az] = distance(lat1, lon1, lat2, lon2);
[latMid, lonMid] = reckon(lat1, lon1, arcDeg/2, az);
end