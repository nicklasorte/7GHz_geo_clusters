function out = buildFinalClustersSimple(cell_data, x_km, n_cluster, cluster_cap)
% buildFinalClustersSimple
%
% Cluster ALL rows into ~n_cluster clusters with a size cap cluster_cap.
% Centers chosen from density ranking (# neighbors within x_km).
% Growth: assign nearest unassigned points to each center until cap.
% Leftovers: assign to nearest center (may overflow cap if infeasible).
%
% Inputs
%   cell_data   : N-by-5 cell {name, lat1, lon1, lat2, lon2} (lat/lon are 1x1 doubles in cells)
%   x_km        : radius for density ranking (used to choose centers)
%   n_cluster   : desired number of clusters
%   cluster_cap : max size per cluster (soft if infeasible)
%
% Output struct out:
%   out.T            : table with name, midpoints, finalClusterID, centerName
%   out.centers      : center indices
%   out.centerNames  : center names
%   out.clusterSizes : sizes after assignment
%   out.info         : diagnostics

T0 = pairsTableFromCell(cell_data);
N = height(T0);

% Midpoints (Mapping Toolbox geodesic midpoint)
[latMid, lonMid] = geodesicMidpoints(T0.lat1, T0.lon1, T0.lat2, T0.lon2);

% Density score for center selection
countWithinKm = countWithinRadius(latMid, lonMid, x_km);
[~, order] = sort(countWithinKm, "descend");

% Choose centers (top of density list)
n_cluster = min(n_cluster, N);
centers = order(1:n_cluster);

% Precompute distances to centers (N-by-n_cluster) in km
D = distancesToCentersKm(latMid, lonMid, centers);

% Greedy fill each cluster up to cap
finalClusterID = zeros(N,1);
centerIndex    = zeros(N,1);
assigned       = false(N,1);

clusterSizes = zeros(n_cluster,1);

for c = 1:n_cluster
    ctr = centers(c);

    % Always include the center itself if unassigned
    if ~assigned(ctr)
        assigned(ctr) = true;
        finalClusterID(ctr) = c;
        centerIndex(ctr) = ctr;
        clusterSizes(c) = clusterSizes(c) + 1;
    end

    % Fill remaining capacity with nearest unassigned points to this center
    space = cluster_cap - clusterSizes(c);
    if space <= 0, continue, end

    d = D(:,c);
    d(assigned) = Inf;
    [~, idxSort] = sort(d, "ascend");

    take = idxSort(1:min(space, sum(isfinite(d))));
    take = take(isfinite(d(take)));

    assigned(take) = true;
    finalClusterID(take) = c;
    centerIndex(take) = ctr;
    clusterSizes(c) = clusterSizes(c) + numel(take);
end

% Assign any leftovers to nearest center (may overflow cap if needed)
left = find(~assigned);
if ~isempty(left)
    [~, nearestC] = min(D(left,:), [], 2);  % nearest center for each leftover
    for k = 1:numel(left)
        i = left(k);
        c = nearestC(k);
        finalClusterID(i) = c;
        centerIndex(i) = centers(c);
        clusterSizes(c) = clusterSizes(c) + 1;
    end
    assigned(left) = true;
end

% Output table (names propagated)
T = T0;
T.latMid = latMid;
T.lonMid = lonMid;
T.countWithinKm = countWithinKm;
T.finalClusterID = finalClusterID;
T.isCenter = false(N,1);
T.isCenter(centers) = true;
T.centerIndex = centerIndex;

T.centerName = strings(N,1);
mask = T.centerIndex > 0;
T.centerName(mask) = T.name(T.centerIndex(mask));

capacity = n_cluster * cluster_cap;
overflow = max(0, N - capacity);

out = struct;
out.T = T;
out.centers = centers;
out.centerNames = T.name(centers);
out.clusterSizes = clusterSizes;
out.info = struct("N",N,"x_km",x_km,"n_cluster",n_cluster,"cluster_cap",cluster_cap, ...
                  "capacity",capacity,"overflowIfInfeasible",overflow);

end

% -------------------------------------------------------------------------
function T = pairsTableFromCell(C)
T = table;
T.name = string(C(:,1));
T.lat1 = cellfun(@(v) v(1), C(:,2));
T.lon1 = cellfun(@(v) v(1), C(:,3));
T.lat2 = cellfun(@(v) v(1), C(:,4));
T.lon2 = cellfun(@(v) v(1), C(:,5));
end

% -------------------------------------------------------------------------
function [latMid, lonMid] = geodesicMidpoints(lat1, lon1, lat2, lon2)
[arcDeg, az] = distance(lat1, lon1, lat2, lon2);
[latMid, lonMid] = reckon(lat1, lon1, arcDeg/2, az);
end

% -------------------------------------------------------------------------
function countWithinKm = countWithinRadius(latMid, lonMid, x_km)
N = numel(latMid);
countWithinKm = zeros(N,1);
for i = 1:N
    dKm = deg2km(distance(latMid(i), lonMid(i), latMid, lonMid));
    countWithinKm(i) = sum(dKm <= x_km) - 1;
end
end

% -------------------------------------------------------------------------
function D = distancesToCentersKm(latMid, lonMid, centers)
% distancesToCentersKm  Compute N-by-C distance matrix (km) to chosen centers.
N = numel(latMid);
C = numel(centers);
D = zeros(N,C);

for c = 1:C
    ctr = centers(c);
    D(:,c) = deg2km(distance(latMid(ctr), lonMid(ctr), latMid, lonMid));
end
end
